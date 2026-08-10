#!/usr/bin/env python3
"""Verify immutable inputs and freeze the primary genetics/state contract."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from genetics_common import (
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    PROJECT_ROOT,
    SCRIPT_DIR,
    assert_candidate_root,
    assert_owned_output,
    atomic_write_tsv,
    bool_text,
    clean,
    ensembl_base,
    inspect_table,
    load_config_tsv,
    parse_float,
    read_table,
    resolve_source,
    sha256_file,
    snapshot_file,
    verify_no_quarantined_inputs,
    write_stage_seal,
)


EXPECTED = {
    "treat_rows": 1918,
    "treat_symbols": 1915,
    "primary_susie_named": 473,
    "primary_susie_bulk_testable": 447,
    "primary_overlap": 34,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    parser.add_argument(
        "--contract",
        type=Path,
        default=SCRIPT_DIR / "config" / "input_contract_v1.tsv",
    )
    return parser.parse_args()


def verify_and_snapshot(
    project: Path, candidate: Path, contract_path: Path
) -> tuple[dict[str, Path], Path]:
    contract = load_config_tsv(contract_path)
    sources: dict[str, Path] = {}
    manifest_rows: list[dict[str, object]] = []

    for item in contract:
        input_id = item["input_id"]
        source = resolve_source(project, item["relative_path"]).resolve()
        if not source.is_file():
            raise ContractError(f"missing input {input_id}: {source}")
        observed_sha = sha256_file(source)
        if observed_sha != item["expected_sha256"]:
            raise ContractError(
                f"hash mismatch for {input_id}: expected {item['expected_sha256']}, "
                f"observed {observed_sha}"
            )

        fmt = item["format"]
        observed_rows = "not_applicable"
        observed_columns = "not_applicable"
        columns_status = "not_applicable"
        expected_rows = int(item["expected_data_rows"])
        if fmt in {"csv", "tsv", "tsv_gz"}:
            header, n_rows = inspect_table(source, fmt)
            observed_rows = str(n_rows)
            observed_columns = ";".join(header)
            required = [x for x in item["required_columns"].split(";") if x]
            absent = sorted(set(required) - set(header))
            if absent:
                raise ContractError(f"{input_id} missing required columns: {absent}")
            columns_status = "pass"
            if expected_rows >= 0 and n_rows != expected_rows:
                raise ContractError(
                    f"row-count mismatch for {input_id}: expected {expected_rows}, observed {n_rows}"
                )

        snapshot_path = ""
        if item["snapshot_mode"] == "copy":
            destination = candidate / "work" / "input_snapshots" / source.name
            assert_owned_output(candidate, destination)
            snapshot_file(source, destination)
            if sha256_file(destination) != observed_sha:
                raise ContractError(f"snapshot hash mismatch for {input_id}")
            snapshot_path = str(destination.relative_to(candidate))
        elif item["snapshot_mode"] != "reference":
            raise ContractError(f"unsupported snapshot mode for {input_id}")

        sources[input_id] = source
        manifest_rows.append(
            {
                "input_id": input_id,
                "source_path": str(source),
                "snapshot_path": snapshot_path,
                "sha256": observed_sha,
                "bytes": source.stat().st_size,
                "expected_data_rows": item["expected_data_rows"],
                "observed_data_rows": observed_rows,
                "required_columns_status": columns_status,
                "observed_columns": observed_columns,
                "snapshot_mode": item["snapshot_mode"],
                "source_role": item["source_role"],
            }
        )

    verify_no_quarantined_inputs(sources.values(), project)
    manifest = candidate / "input_manifest.tsv"
    assert_owned_output(candidate, manifest)
    fields = [
        "input_id",
        "source_path",
        "snapshot_path",
        "sha256",
        "bytes",
        "expected_data_rows",
        "observed_data_rows",
        "required_columns_status",
        "observed_columns",
        "snapshot_mode",
        "source_role",
    ]
    atomic_write_tsv(manifest, manifest_rows, fields)
    return sources, manifest


def collapse_bulk(rows: list[dict[str, str]]) -> tuple[dict[str, dict[str, object]], int]:
    by_symbol: dict[str, list[dict[str, str]]] = defaultdict(list)
    treat_rows = 0
    for row in rows:
        symbol = clean(row.get("symbol"))
        treat_fdr = parse_float(row.get("treat_fdr"))
        if treat_fdr is not None and treat_fdr < 0.05:
            treat_rows += 1
        if symbol:
            by_symbol[symbol].append(row)

    collapsed: dict[str, dict[str, object]] = {}
    for symbol, candidates in by_symbol.items():
        ranked = sorted(
            candidates,
            key=lambda row: (
                parse_float(row.get("treat_fdr"))
                if parse_float(row.get("treat_fdr")) is not None
                else float("inf"),
                -(abs(parse_float(row.get("t")) or 0.0)),
                clean(row.get("gene")),
            ),
        )
        selected = ranked[0]
        collapsed[symbol] = {
            "ensembl_bulk": ";".join(
                sorted({ensembl_base(row.get("gene")) for row in candidates if ensembl_base(row.get("gene"))})
            ),
            "bulk_row_count": len(candidates),
            "bulk_logFC": clean(selected.get("logFC")),
            "bulk_t": clean(selected.get("t")),
            "bulk_AveExpr": clean(selected.get("AveExpr")),
            "bulk_treat_fdr": clean(selected.get("treat_fdr")),
            "established_state_associated": any(
                parse_float(row.get("treat_fdr")) is not None
                and float(row["treat_fdr"]) < 0.05
                for row in candidates
            ),
        }
    return collapsed, treat_rows


def collapse_genetics(
    tier_rows: list[dict[str, str]], full_rows: list[dict[str, str]]
) -> dict[str, dict[str, object]]:
    ensembl_to_symbols: dict[str, set[str]] = defaultdict(set)
    for row in full_rows:
        ens = ensembl_base(row.get("ensembl"))
        symbol = clean(row.get("gene"))
        if ens and symbol:
            ensembl_to_symbols[ens].add(symbol)

    by_symbol: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tier_rows:
        ens = ensembl_base(row.get("ensembl"))
        if not ens:
            continue
        symbols = ensembl_to_symbols.get(ens, set())
        if len(symbols) > 1:
            raise ContractError(f"ambiguous full-gene symbol map for {ens}: {sorted(symbols)}")
        if not symbols:
            continue
        enriched = dict(row)
        enriched["ensembl_base"] = ens
        by_symbol[next(iter(symbols))].append(enriched)

    collapsed: dict[str, dict[str, object]] = {}
    for symbol, candidates in by_symbol.items():
        ranked = sorted(
            candidates,
            key=lambda row: (
                -(
                    parse_float(row.get("coloc_best_susie_pp4"))
                    if parse_float(row.get("coloc_best_susie_pp4")) is not None
                    else -1.0
                ),
                row["ensembl_base"],
            ),
        )
        selected = ranked[0]
        max_susie = max(
            (parse_float(row.get("coloc_best_susie_pp4")) for row in candidates),
            default=None,
            key=lambda value: -1.0 if value is None else value,
        )
        max_abf = max(
            (parse_float(row.get("coloc_best_abf_pp4")) for row in candidates),
            default=None,
            key=lambda value: -1.0 if value is None else value,
        )
        collapsed[symbol] = {
            "ensembl_genetic": ";".join(sorted({row["ensembl_base"] for row in candidates})),
            "genetic_row_count": len(candidates),
            "coloc_best_susie_pp4": "" if max_susie is None else max_susie,
            "coloc_best_abf_pp4": "" if max_abf is None else max_abf,
            "driving_gwas": clean(selected.get("driving_gwas")),
            "driving_trait": clean(selected.get("driving_trait")),
            "primary_genetic": max_susie is not None and max_susie > 0.5,
        }
    return collapsed


def build_frozen_classes(
    sources: dict[str, Path], candidate: Path
) -> tuple[Path, Path]:
    bulk_rows = read_table(sources["canonical_treat"], "csv")
    tier_rows = read_table(sources["tier12_genelevel"], "csv")
    full_rows = read_table(sources["full_genelevel"], "csv")
    bulk, treat_rows = collapse_bulk(bulk_rows)
    genetics = collapse_genetics(tier_rows, full_rows)

    all_symbols = sorted(set(bulk) | set(genetics))
    frozen_rows: list[dict[str, object]] = []
    for symbol in all_symbols:
        b = bulk.get(symbol, {})
        g = genetics.get(symbol, {})
        bulk_tested = symbol in bulk
        genetic_tested = symbol in genetics
        joint = bulk_tested and genetic_tested
        primary = bool(g.get("primary_genetic", False))
        state = bool(b.get("established_state_associated", False))
        if not joint:
            static_class = "indeterminate_not_jointly_testable"
        elif primary and state:
            static_class = "convergent"
        elif primary:
            static_class = "genetic_only"
        elif state:
            static_class = "disease_state_only"
        else:
            static_class = "neither"
        frozen_rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_bulk": b.get("ensembl_bulk", ""),
                "ensembl_genetic": g.get("ensembl_genetic", ""),
                "bulk_tested": bool_text(bulk_tested),
                "primary_genetic_map_tested": bool_text(genetic_tested),
                "joint_testable": bool_text(joint),
                "primary_genetic": bool_text(primary),
                "established_state_associated": bool_text(state),
                "static_class": static_class,
                "bulk_row_count": b.get("bulk_row_count", 0),
                "bulk_logFC": b.get("bulk_logFC", ""),
                "bulk_t": b.get("bulk_t", ""),
                "bulk_AveExpr": b.get("bulk_AveExpr", ""),
                "bulk_treat_fdr": b.get("bulk_treat_fdr", ""),
                "genetic_row_count": g.get("genetic_row_count", 0),
                "coloc_best_susie_pp4": g.get("coloc_best_susie_pp4", ""),
                "coloc_best_abf_pp4": g.get("coloc_best_abf_pp4", ""),
                "driving_gwas": g.get("driving_gwas", ""),
                "driving_trait": g.get("driving_trait", ""),
                "context_annotation": "not_evaluated_preflight",
            }
        )

    observed = {
        "treat_rows": treat_rows,
        "treat_symbols": sum(
            1 for value in bulk.values() if value["established_state_associated"]
        ),
        "primary_susie_named": sum(1 for value in genetics.values() if value["primary_genetic"]),
        "primary_susie_bulk_testable": sum(
            1
            for symbol, value in genetics.items()
            if value["primary_genetic"] and symbol in bulk
        ),
        "primary_overlap": sum(
            1
            for symbol, value in genetics.items()
            if value["primary_genetic"]
            and symbol in bulk
            and bulk[symbol]["established_state_associated"]
        ),
    }
    audit_rows = []
    failures = []
    for metric, expected in EXPECTED.items():
        actual = observed[metric]
        passed = actual == expected
        audit_rows.append(
            {
                "metric": metric,
                "expected": expected,
                "observed": actual,
                "status": "pass" if passed else "fail",
                "definition": {
                    "treat_rows": "canonical rows with treat_fdr < 0.05",
                    "treat_symbols": "unique nonempty symbols with any treat_fdr < 0.05",
                    "primary_susie_named": "unique named genes with SuSiE PP.H4 > 0.5",
                    "primary_susie_bulk_testable": "primary genes with a nonempty canonical bulk symbol",
                    "primary_overlap": "primary bulk-testable genes with TREAT FDR < 0.05",
                }[metric],
            }
        )
        if not passed:
            failures.append(f"{metric}: expected {expected}, observed {actual}")
    if failures:
        raise ContractError("frozen contract failed: " + "; ".join(failures))

    classes_path = candidate / "frozen_evidence_classes.tsv"
    audit_path = candidate / "frozen_contract_rederivation.tsv"
    fields = [
        "gene_symbol",
        "ensembl_bulk",
        "ensembl_genetic",
        "bulk_tested",
        "primary_genetic_map_tested",
        "joint_testable",
        "primary_genetic",
        "established_state_associated",
        "static_class",
        "bulk_row_count",
        "bulk_logFC",
        "bulk_t",
        "bulk_AveExpr",
        "bulk_treat_fdr",
        "genetic_row_count",
        "coloc_best_susie_pp4",
        "coloc_best_abf_pp4",
        "driving_gwas",
        "driving_trait",
        "context_annotation",
    ]
    atomic_write_tsv(classes_path, frozen_rows, fields)
    atomic_write_tsv(
        audit_path,
        audit_rows,
        ["metric", "expected", "observed", "status", "definition"],
    )
    return classes_path, audit_path


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    candidate = assert_candidate_root(project, args.candidate_root)
    candidate.mkdir(parents=True, exist_ok=True)
    sources, manifest = verify_and_snapshot(project, candidate, args.contract.resolve())
    classes, audit = build_frozen_classes(sources, candidate)
    write_stage_seal(
        candidate,
        "01_freeze_and_rederive",
        [manifest, classes, audit],
        [args.contract.resolve()],
    )
    print(f"PASS: frozen contract written beneath {candidate}")


if __name__ == "__main__":
    main()
