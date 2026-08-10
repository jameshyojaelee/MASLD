#!/usr/bin/env python3
"""Outcome-blind Ensembl-v49 identity adjudication for frozen GEN classes.

The producer never chooses among multiple valid Ensembl IDs using an effect,
P value, testability flag, or evidence class.  It emits one accepted row per
stable GENCODE-v49 ID, preserves every contributing source row verbatim as
canonical JSON, and quarantines identity ambiguity or conflicting duplicate
calls.  The products are candidate-only inputs for PASS; they are not a
canonical GEN release and cannot authorize manuscript promotion.
"""

from __future__ import annotations

import argparse
import json
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd

import generate_evidence_passports as contract


SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_OUTPUT_ROOT = (
    contract.DEFAULT_CANDIDATE_ROOT / "preflight" / "gen_identity_adjudication_v1"
)
ENSG_PATTERN = re.compile(r"^ENSG\d{11}$")
IDENTITY_COLUMNS = {"gene_symbol", "ensembl_bulk", "ensembl_genetic"}
ADJUDICATION_COLUMNS = [
    "adjudicated_ensembl_id",
    "adjudicated_symbol",
    "identity_resolution_rule",
    "source_row_count",
    "source_row_ids",
    "source_row_sha256s",
    "source_rows_json",
    "source_call_sha256",
    "identity_adjudication_status",
]
QUARANTINE_COLUMNS = [
    "quarantine_id",
    "source_row_id",
    "source_row_sha256",
    "source_symbol",
    "source_ensembl_bulk",
    "source_ensembl_genetic",
    "normalized_ensembl_candidates",
    "reason_code",
    "reason",
    "source_artifact_sha256",
    "source_row_json",
]
AUDIT_COLUMNS = [
    "release_id",
    "status",
    "source_artifact_path",
    "source_artifact_sha256",
    "identity_artifact_path",
    "identity_artifact_sha256",
    "n_source_rows",
    "n_individually_resolved_rows",
    "n_adjudicated_genes",
    "n_collapsed_duplicate_rows",
    "n_quarantined_source_rows",
    "n_irreducible_multi_id",
    "n_absent_v49",
    "n_identity_conflict",
    "n_unresolved",
    "n_invalid_token",
    "outcome_fields_used",
    "merge_discordant_calls",
    "one_row_per_ensembl",
    "adjudicated_artifact_sha256",
    "quarantine_artifact_sha256",
    "producer_sha256",
]
READY_COLUMNS = [
    "release_id",
    "status",
    "source_artifact_sha256",
    "identity_artifact_sha256",
    "adjudicated_artifact_sha256",
    "quarantine_artifact_sha256",
    "audit_artifact_sha256",
    "producer_sha256",
    "n_source_rows",
    "n_adjudicated_genes",
    "n_quarantined_source_rows",
    "n_irreducible_multi_id",
    "n_absent_v49",
    "n_identity_conflict",
    "n_unresolved",
    "n_invalid_token",
    "outcome_fields_used",
    "merge_discordant_calls",
    "one_row_per_ensembl",
    "canonical_promotion_authorized",
    "fixture_only",
]


def _clean(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path, engine="pyarrow").astype("string").fillna("")
    compression = "gzip" if path.suffix == ".gz" else "infer"
    return pd.read_csv(
        path, sep="\t", dtype="string", keep_default_na=False,
        compression=compression,
    )


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _row_hash(row: Mapping[str, Any]) -> str:
    return contract.stable_sha256(
        {"row": {str(key): _clean(value) for key, value in sorted(row.items())}}
    )


def _write_tsv(path: Path, rows: Iterable[Mapping[str, Any]], columns: list[str]) -> None:
    contract.atomic_write_tsv(path, rows, columns)


def _load_identity(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    frame = _read_table(path)
    id_column = "ensembl_base" if "ensembl_base" in frame else "gene_id"
    symbol_column = "gene_name" if "gene_name" in frame else "symbol"
    missing = {id_column, symbol_column} - set(frame.columns)
    if missing:
        raise contract.PassportContractError(
            "GEN_IDENTITY_SCHEMA", f"missing GENCODE columns: {sorted(missing)}"
        )
    by_id: dict[str, str] = {}
    for record in frame.to_dict("records"):
        stable_id = _clean(record[id_column]).split(".", 1)[0]
        symbol = _clean(record[symbol_column])
        if not ENSG_PATTERN.fullmatch(stable_id) or not symbol:
            continue
        previous = by_id.get(stable_id)
        if previous is not None and previous != symbol:
            raise contract.PassportContractError(
                "GEN_IDENTITY_REFERENCE_CONFLICT",
                f"{stable_id}: {previous!r} versus {symbol!r}",
            )
        by_id[stable_id] = symbol
    by_symbol: dict[str, list[str]] = defaultdict(list)
    for stable_id, symbol in sorted(by_id.items()):
        by_symbol[symbol].append(stable_id)
    if not by_id:
        raise contract.PassportContractError("GEN_IDENTITY_REFERENCE_EMPTY", str(path))
    return by_id, dict(by_symbol)


def _parse_identity_tokens(values: Iterable[Any]) -> tuple[set[str], list[str]]:
    candidates: set[str] = set()
    invalid: list[str] = []
    for value in values:
        raw = _clean(value)
        if not raw:
            continue
        for token in re.split(r"[;,|]", raw):
            normalized = token.strip().split(".", 1)[0]
            if not normalized:
                continue
            if ENSG_PATTERN.fullmatch(normalized):
                candidates.add(normalized)
            else:
                invalid.append(normalized)
    return candidates, sorted(set(invalid))


def _quarantine_row(
    source: Mapping[str, Any],
    row_id: str,
    row_sha: str,
    candidates: Iterable[str],
    code: str,
    reason: str,
    source_sha: str,
) -> dict[str, str]:
    return {
        "quarantine_id": "gen_identity:" + contract.stable_sha256(
            {"row_id": row_id, "reason_code": code, "source_row_sha256": row_sha}
        )[:20],
        "source_row_id": row_id,
        "source_row_sha256": row_sha,
        "source_symbol": _clean(source.get("gene_symbol")),
        "source_ensembl_bulk": _clean(source.get("ensembl_bulk")),
        "source_ensembl_genetic": _clean(source.get("ensembl_genetic")),
        "normalized_ensembl_candidates": ";".join(sorted(set(candidates))),
        "reason_code": code,
        "reason": reason,
        "source_artifact_sha256": source_sha,
        "source_row_json": _canonical_json(
            {str(key): _clean(value) for key, value in source.items()}
        ),
    }


def adjudicate(
    source_path: Path,
    identity_path: Path,
    output_root: Path,
    release_id: str,
    fixture_only: bool,
) -> dict[str, Path]:
    """Create deterministic adjudicated, quarantine, audit, and READY products."""

    source_path = source_path.resolve()
    identity_path = identity_path.resolve()
    output_root = output_root.resolve()
    if not fixture_only and output_root != DEFAULT_OUTPUT_ROOT.resolve():
        raise contract.PassportContractError(
            "GEN_ADJUDICATION_OUTPUT_ROOT",
            f"production candidate output must be {DEFAULT_OUTPUT_ROOT.resolve()}",
        )
    if not source_path.is_file() or not identity_path.is_file():
        raise contract.PassportContractError(
            "GEN_ADJUDICATION_INPUT_MISSING", f"source={source_path}; identity={identity_path}"
        )
    source = _read_table(source_path)
    required = {
        "gene_symbol", "ensembl_bulk", "ensembl_genetic", "static_class",
        "bulk_tested", "primary_genetic_map_tested", "primary_genetic",
        "established_state_associated",
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise contract.PassportContractError(
            "GEN_ADJUDICATION_SOURCE_SCHEMA", f"missing={missing}"
        )
    by_id, by_symbol = _load_identity(identity_path)
    valid_identity_ids = set(by_id)
    source_sha = contract.sha256_file(source_path)
    identity_sha = contract.sha256_file(identity_path)

    resolved: dict[str, list[dict[str, Any]]] = defaultdict(list)
    quarantine: list[dict[str, str]] = []
    reason_counts: Counter[str] = Counter()
    individually_resolved = 0
    for index, raw_record in enumerate(source.to_dict("records"), start=1):
        record = {str(key): _clean(value) for key, value in raw_record.items()}
        row_id = f"GEN_SOURCE:{index:06d}"
        row_sha = _row_hash(record)
        candidates, invalid = _parse_identity_tokens(
            [record.get("ensembl_bulk", ""), record.get("ensembl_genetic", "")]
        )
        valid = sorted(candidates & valid_identity_ids)
        rule = ""
        stable_id = ""
        if invalid:
            code = "GEN_INVALID_ENSEMBL_TOKEN"
            reason = "At least one deposited identity token is not a stable human Ensembl gene ID."
        elif len(valid) > 1:
            code = "GEN_IRREDUCIBLE_MULTI_ENSEMBL"
            reason = "More than one deposited candidate is a valid GENCODE-v49 stable ID; no outcome-blind identity choice is possible."
        elif len(valid) == 1:
            stable_id = valid[0]
            rule = "unique_deposited_v49_candidate"
            code = ""
            reason = ""
        elif candidates:
            code = "GEN_ENSEMBL_ABSENT_V49"
            reason = "Deposited stable Ensembl candidate(s) are absent from the selected GENCODE-v49 identity reference."
        else:
            symbol_matches = by_symbol.get(record.get("gene_symbol", ""), [])
            if len(symbol_matches) == 1:
                stable_id = symbol_matches[0]
                rule = "unique_v49_symbol_fallback_no_deposited_ensembl"
                code = ""
                reason = ""
            elif len(symbol_matches) > 1:
                code = "GEN_IRREDUCIBLE_MULTI_ENSEMBL"
                reason = "No Ensembl ID was deposited and the source symbol maps to multiple GENCODE-v49 stable IDs."
            else:
                code = "GEN_ENSEMBL_UNRESOLVED"
                reason = "Neither a deposited Ensembl ID nor a unique GENCODE-v49 symbol mapping is available."
        if code:
            quarantine.append(
                _quarantine_row(
                    record, row_id, row_sha, candidates, code, reason, source_sha
                )
            )
            reason_counts[code] += 1
            continue
        individually_resolved += 1
        resolved[stable_id].append(
            {
                "record": record,
                "row_id": row_id,
                "row_sha": row_sha,
                "resolution_rule": rule,
            }
        )

    adjudicated: list[dict[str, Any]] = []
    collapsed_duplicates = 0
    for stable_id in sorted(resolved):
        entries = resolved[stable_id]
        call_payloads = [
            {
                key: value
                for key, value in entry["record"].items()
                if key not in IDENTITY_COLUMNS
            }
            for entry in entries
        ]
        call_hashes = {_row_hash(payload) for payload in call_payloads}
        if len(call_hashes) != 1:
            reason = (
                "Multiple source rows resolve to the same stable Ensembl ID but carry "
                "non-identical source calls; PASS will not choose or merge them."
            )
            for entry in entries:
                quarantine.append(
                    _quarantine_row(
                        entry["record"], entry["row_id"], entry["row_sha"],
                        [stable_id], "GEN_IDENTITY_CONFLICT", reason, source_sha,
                    )
                )
                reason_counts["GEN_IDENTITY_CONFLICT"] += 1
            continue
        first = entries[0]
        source_rows = [entry["record"] for entry in entries]
        row = dict(first["record"])
        row.update(
            {
                "adjudicated_ensembl_id": stable_id,
                "adjudicated_symbol": by_id[stable_id],
                "identity_resolution_rule": ";".join(
                    sorted({entry["resolution_rule"] for entry in entries})
                ),
                "source_row_count": str(len(entries)),
                "source_row_ids": ";".join(entry["row_id"] for entry in entries),
                "source_row_sha256s": ";".join(entry["row_sha"] for entry in entries),
                "source_rows_json": _canonical_json(source_rows),
                "source_call_sha256": next(iter(call_hashes)),
                "identity_adjudication_status": "resolved_v49",
            }
        )
        adjudicated.append(row)
        collapsed_duplicates += len(entries) - 1

    if len({row["adjudicated_ensembl_id"] for row in adjudicated}) != len(adjudicated):
        raise contract.PassportContractError(
            "GEN_ADJUDICATION_DUPLICATE_OUTPUT", "adjudicated output is not one row per stable ID"
        )
    output_root.mkdir(parents=True, exist_ok=True)
    adjudicated_path = output_root / "gen_frozen_classes_ensembl_adjudicated.tsv"
    quarantine_path = output_root / "gen_identity_quarantine.tsv"
    audit_path = output_root / "gen_identity_adjudication_audit.tsv"
    ready_path = output_root / "GEN_IDENTITY_ADJUDICATION_READY"
    _write_tsv(
        adjudicated_path, adjudicated, [*list(source.columns), *ADJUDICATION_COLUMNS]
    )
    _write_tsv(quarantine_path, quarantine, QUARANTINE_COLUMNS)
    n_source = len(source)
    n_quarantine = len(quarantine)
    audit_row = {
        "release_id": release_id,
        "status": "identity_adjudication_complete_with_quarantine",
        "source_artifact_path": str(source_path),
        "source_artifact_sha256": source_sha,
        "identity_artifact_path": str(identity_path),
        "identity_artifact_sha256": identity_sha,
        "n_source_rows": n_source,
        "n_individually_resolved_rows": individually_resolved,
        "n_adjudicated_genes": len(adjudicated),
        "n_collapsed_duplicate_rows": collapsed_duplicates,
        "n_quarantined_source_rows": n_quarantine,
        "n_irreducible_multi_id": reason_counts["GEN_IRREDUCIBLE_MULTI_ENSEMBL"],
        "n_absent_v49": reason_counts["GEN_ENSEMBL_ABSENT_V49"],
        "n_identity_conflict": reason_counts["GEN_IDENTITY_CONFLICT"],
        "n_unresolved": reason_counts["GEN_ENSEMBL_UNRESOLVED"],
        "n_invalid_token": reason_counts["GEN_INVALID_ENSEMBL_TOKEN"],
        "outcome_fields_used": "false",
        "merge_discordant_calls": "false",
        "one_row_per_ensembl": "true",
        "adjudicated_artifact_sha256": contract.sha256_file(adjudicated_path),
        "quarantine_artifact_sha256": contract.sha256_file(quarantine_path),
        "producer_sha256": contract.sha256_file(SCRIPT_PATH),
    }
    _write_tsv(audit_path, [audit_row], AUDIT_COLUMNS)
    ready_row = {
        "release_id": release_id,
        "status": "identity_adjudication_complete_with_quarantine",
        "source_artifact_sha256": source_sha,
        "identity_artifact_sha256": identity_sha,
        "adjudicated_artifact_sha256": contract.sha256_file(adjudicated_path),
        "quarantine_artifact_sha256": contract.sha256_file(quarantine_path),
        "audit_artifact_sha256": contract.sha256_file(audit_path),
        "producer_sha256": contract.sha256_file(SCRIPT_PATH),
        "n_source_rows": n_source,
        "n_adjudicated_genes": len(adjudicated),
        "n_quarantined_source_rows": n_quarantine,
        "n_irreducible_multi_id": reason_counts["GEN_IRREDUCIBLE_MULTI_ENSEMBL"],
        "n_absent_v49": reason_counts["GEN_ENSEMBL_ABSENT_V49"],
        "n_identity_conflict": reason_counts["GEN_IDENTITY_CONFLICT"],
        "n_unresolved": reason_counts["GEN_ENSEMBL_UNRESOLVED"],
        "n_invalid_token": reason_counts["GEN_INVALID_ENSEMBL_TOKEN"],
        "outcome_fields_used": "false",
        "merge_discordant_calls": "false",
        "one_row_per_ensembl": "true",
        "canonical_promotion_authorized": "false",
        "fixture_only": str(fixture_only).lower(),
    }
    _write_tsv(ready_path, [ready_row], READY_COLUMNS)
    return {
        "adjudicated": adjudicated_path,
        "quarantine": quarantine_path,
        "audit": audit_path,
        "ready": ready_path,
    }


def run_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="passport_gen_identity_") as temp_dir:
        root = Path(temp_dir)
        identity = root / "identity.tsv"
        source = root / "source.tsv"
        _write_tsv(
            identity,
            [
                {"ensembl_base": "ENSG00000000001", "gene_name": "A"},
                {"ensembl_base": "ENSG00000000002", "gene_name": "B"},
                {"ensembl_base": "ENSG00000000003", "gene_name": "C"},
            ],
            ["ensembl_base", "gene_name"],
        )
        base = {
            "bulk_tested": "true", "primary_genetic_map_tested": "true",
            "primary_genetic": "false", "established_state_associated": "true",
            "static_class": "disease_state_only", "bulk_logFC": "0.5",
        }
        rows = [
            {"gene_symbol": "A", "ensembl_bulk": "ENSG00000000001", "ensembl_genetic": "", **base},
            {"gene_symbol": "B", "ensembl_bulk": "ENSG00000000002;ENSG00000000003", "ensembl_genetic": "", **base},
            {"gene_symbol": "OLD", "ensembl_bulk": "ENSG99999999999", "ensembl_genetic": "", **base},
            {"gene_symbol": "C", "ensembl_bulk": "ENSG00000000003", "ensembl_genetic": "", **base},
            {"gene_symbol": "C", "ensembl_bulk": "ENSG00000000003", "ensembl_genetic": "", **{**base, "bulk_logFC": "0.8"}},
        ]
        _write_tsv(source, rows, list(rows[0]))
        products = adjudicate(source, identity, root / "out", "fixture", True)
        accepted = contract.read_tsv_rows(products["adjudicated"])
        quarantined = contract.read_tsv_rows(products["quarantine"])
        ready = contract.read_tsv_rows(products["ready"])[0]
        if [row["adjudicated_ensembl_id"] for row in accepted] != ["ENSG00000000001"]:
            raise contract.PassportContractError(
                "GEN_ADJUDICATION_SELF_TEST", "unexpected accepted stable-ID set"
            )
        observed = Counter(row["reason_code"] for row in quarantined)
        expected = Counter(
            {
                "GEN_IRREDUCIBLE_MULTI_ENSEMBL": 1,
                "GEN_ENSEMBL_ABSENT_V49": 1,
                "GEN_IDENTITY_CONFLICT": 2,
            }
        )
        if observed != expected or ready["outcome_fields_used"] != "false":
            raise contract.PassportContractError(
                "GEN_ADJUDICATION_SELF_TEST", f"observed={dict(observed)}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--identity", type=Path)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--release-id", default=contract.PRODUCTION_ANALYSIS_RELEASE_ID
    )
    parser.add_argument("--fixture-mode", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_test()
        print("PASS: GEN Ensembl identity adjudication self-test")
        return
    if args.source is None or args.identity is None:
        raise contract.PassportContractError(
            "GEN_ADJUDICATION_ARGUMENT", "--source and --identity are required"
        )
    products = adjudicate(
        args.source, args.identity, args.output_root, args.release_id, args.fixture_mode
    )
    print(
        "PASS: candidate-only GEN identity adjudication; "
        + ";".join(f"{name}={path}" for name, path in products.items())
    )


if __name__ == "__main__":
    try:
        main()
    except contract.PassportContractError as exc:
        raise SystemExit(str(exc)) from exc
