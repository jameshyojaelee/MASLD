#!/usr/bin/env python3
"""Verify immutable inputs and freeze the primary genetics/state contract."""

from __future__ import annotations

import argparse
import json
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
    iter_table,
    load_config_tsv,
    parse_float,
    read_table,
    resolve_source,
    sha256_file,
    snapshot_file,
    verify_no_quarantined_inputs,
    write_stage_seal,
)


# Method-specific ABF driver columns written by rebuild_tier12.py (review item 5).
# They are carried only when the Tier-1/2 table has them.
ABF_DRIVER_FIELDS = ("abf_driving_gwas", "abf_driving_trait")

# Counts recorded by the v1 freeze (2026-08-07). They are a historical record and
# the template of the metric names; they are not asserted, because the shared
# bulk-row rule and the Ensembl join (2026-09-24) move them on the v1 inputs
# (treat_symbols 1,909, primary_susie_bulk_testable 449, primary_overlap 33).
# Every run passes --expected, computed from the inputs by
# 11_write_corrected_input_contract.py.
EXPECTED = {
    "treat_rows": 1918,
    "treat_symbols": 1915,
    "primary_susie_named": 473,
    "primary_susie_bulk_testable": 447,
    "primary_overlap": 34,
}

# Established-state rules. "treat" (default) is the rule of the frozen 2026-08-07
# contract. "canonical" is the paper's DEG rule since 2026-08-12 and adds
# bulk_padj and established_state_rule columns.
STATE_RULES = {
    "treat": "treat_fdr0.05",
    "canonical": "canonical_padj0.05_absLFC0.5",
}
# Corrected 06_susie_coloc.R: every SuSiE signal pair failed coloc's
# shared-posterior check. The gene is untestable by SuSiE, never a negative.
UNTESTABLE_METHOD = "susie_untestable_insufficient_shared_posterior"
UNTESTABLE_CLASS = "genetic_untestable_shared_posterior"
# Shared bulk-row rule, identical to collapse_bulk_rows() in
# scripts/manuscript/build_evidence_class_release.R: when several bulk rows carry
# one symbol, the row on the primary assembly (GENCODE v49 chr1-22, X, Y, M) is
# kept; alternate-haplotype and patch rows (e.g. HLA-DRA on GL000252.2) are used
# only when the symbol has no primary row. Ties within an assembly class go to
# the lowest q-value of the state rule, then the larger |t|, then the gene ID.
# The kept row alone gives the state call and the displayed values.
PRIMARY_CHROMOSOMES = frozenset([f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY", "chrM"])


def state_positive(row: dict[str, str], rule: str) -> bool:
    if rule == "treat":
        treat_fdr = parse_float(row.get("treat_fdr"))
        return treat_fdr is not None and treat_fdr < 0.05
    padj = parse_float(row.get("padj"))
    logfc = parse_float(row.get("logFC"))
    return padj is not None and logfc is not None and padj < 0.05 and abs(logfc) > 0.5


def metric_names(rule: str, untestable: bool) -> list[str]:
    """Frozen-contract count names; the v1 names when rule=treat without untestable."""
    names = [f"{rule}_rows", f"{rule}_symbols", *list(EXPECTED)[2:]]
    return names + (["genetic_untestable_named"] if untestable else [])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    parser.add_argument(
        "--contract",
        type=Path,
        default=SCRIPT_DIR / "config" / "input_contract_v1.tsv",
    )
    parser.add_argument(
        "--expected",
        type=Path,
        required=True,
        help="JSON of the frozen-contract counts (11_write_corrected_input_contract.py)",
    )
    parser.add_argument("--state-rule", choices=sorted(STATE_RULES), default="treat")
    parser.add_argument(
        "--untestable-state",
        action="store_true",
        help=f"give genes whose only Tier-1/2 SuSiE result is {UNTESTABLE_METHOD} the class {UNTESTABLE_CLASS}",
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


def collapse_bulk(
    rows: list[dict[str, str]], rule: str = "treat", primary_ids: set[str] | None = None
) -> tuple[dict[str, dict[str, object]], int]:
    by_symbol: dict[str, list[dict[str, str]]] = defaultdict(list)
    state_rows = 0
    # The displayed row is the most significant one under the rule's own q-value.
    q_field = "treat_fdr" if rule == "treat" else "padj"
    for row in rows:
        symbol = clean(row.get("symbol"))
        if state_positive(row, rule):
            state_rows += 1
        if symbol:
            by_symbol[symbol].append(row)

    collapsed: dict[str, dict[str, object]] = {}
    for symbol, candidates in by_symbol.items():
        selected = select_bulk_row(candidates, q_field, primary_ids)
        # Identity candidates: the primary-assembly IDs when the symbol has any
        # (outcome-blind; alternate-contig copies drop out), else every ID.
        ids = {ensembl_base(row.get("gene")) for row in candidates if ensembl_base(row.get("gene"))}
        if primary_ids is not None and ids & primary_ids:
            ids &= primary_ids
        collapsed[symbol] = {
            "ensembl_bulk": ";".join(sorted(ids)),
            "bulk_row_count": len(candidates),
            "bulk_logFC": clean(selected.get("logFC")),
            "bulk_t": clean(selected.get("t")),
            "bulk_AveExpr": clean(selected.get("AveExpr")),
            "bulk_treat_fdr": clean(selected.get("treat_fdr")),
            "bulk_padj": clean(selected.get("padj")),
            "bulk_selected_ensembl": ensembl_base(selected.get("gene")),
            "established_state_associated": state_positive(selected, rule),
        }
    return collapsed, state_rows


def primary_assembly_ids(metadata_path: Path) -> set[str]:
    """Unversioned Ensembl IDs on the primary assembly in GENCODE v49 metadata."""
    return {
        ensembl_base(row.get("gene_id"))
        for row in iter_table(metadata_path, "tsv_gz")
        if clean(row.get("chromosome")) in PRIMARY_CHROMOSOMES
    }


def select_bulk_row(
    candidates: list[dict[str, str]], q_field: str, primary_ids: set[str] | None = None
) -> dict[str, str]:
    """Shared bulk-row rule; see PRIMARY_CHROMOSOMES."""
    return min(
        candidates,
        key=lambda row: (
            primary_ids is not None and ensembl_base(row.get("gene")) not in primary_ids,
            parse_float(row.get(q_field)) if parse_float(row.get(q_field)) is not None else float("inf"),
            -(abs(parse_float(row.get("t")) or 0.0)),
            clean(row.get("gene")),
        ),
    )


def untestable_ensembl(all_gwas_path: Path, tier_path: Path) -> set[str]:
    """Ensembl IDs with at least one Tier-1/2 pair of the untestable SuSiE method."""
    tier12 = {
        row["study_name"] for row in read_table(tier_path, "tsv") if clean(row.get("tier")) in {"1", "2"}
    }
    found: set[str] = set()
    for row in iter_table(all_gwas_path, "csv"):
        if row.get("method") == UNTESTABLE_METHOD and clean(row.get("gwas_name")) in tier12:
            ens = ensembl_base(row.get("ensembl"))
            if ens:
                found.add(ens)
    return found


def collapse_genetics(
    tier_rows: list[dict[str, str]],
    full_rows: list[dict[str, str]],
    untestable: set[str] | None = None,
    bulk_symbols: dict[str, str] | None = None,
) -> dict[str, dict[str, object]]:
    """Collapse Tier-1/2 COLOC rows to symbols.

    A COLOC gene joins the bulk table by unversioned Ensembl ID: when its ID
    carries a bulk row, the bulk (GENCODE v49) symbol is used, so renamed genes
    (C11orf80 -> TOP6BL) stay one gene. Other IDs keep the COLOC symbol.
    """
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
        if bulk_symbols is not None and ens in bulk_symbols:
            symbols = {bulk_symbols[ens]}
        else:
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
        if untestable is not None:
            collapsed[symbol]["genetic_untestable"] = not collapsed[symbol]["primary_genetic"] and any(
                row["ensembl_base"] in untestable for row in candidates
            )
        if all(field in selected for field in ABF_DRIVER_FIELDS):
            # The ABF driver comes from the row holding the maximum ABF PP4.
            abf_selected = sorted(
                candidates,
                key=lambda row: (
                    -(
                        parse_float(row.get("coloc_best_abf_pp4"))
                        if parse_float(row.get("coloc_best_abf_pp4")) is not None
                        else -1.0
                    ),
                    row["ensembl_base"],
                ),
            )[0]
            for field in ABF_DRIVER_FIELDS:
                collapsed[symbol][field] = clean(abf_selected.get(field))
    return collapsed


def build_frozen_classes(
    sources: dict[str, Path],
    candidate: Path,
    expected_counts: dict[str, int] = EXPECTED,
    rule: str = "treat",
    with_untestable: bool = False,
) -> tuple[Path, Path]:
    bulk_rows = read_table(sources["canonical_treat"], "csv")
    tier_rows = read_table(sources["tier12_genelevel"], "csv")
    full_rows = read_table(sources["full_genelevel"], "csv")
    bulk, state_rows = collapse_bulk(
        bulk_rows, rule, primary_assembly_ids(sources["gencode_metadata"])
    )
    bulk_symbols = {
        ensembl_base(row.get("gene")): clean(row.get("symbol"))
        for row in bulk_rows
        if ensembl_base(row.get("gene")) and clean(row.get("symbol"))
    }
    untestable = (
        untestable_ensembl(sources["all_gwas_coloc"], sources["gwas_trait_tier"])
        if with_untestable
        else None
    )
    genetics = collapse_genetics(tier_rows, full_rows, untestable, bulk_symbols)
    has_abf_driver = bool(tier_rows) and all(field in tier_rows[0] for field in ABF_DRIVER_FIELDS)

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
        if g.get("genetic_untestable", False):
            static_class = UNTESTABLE_CLASS
        elif not joint:
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
        if has_abf_driver:
            frozen_rows[-1].update({field: g.get(field, "") for field in ABF_DRIVER_FIELDS})
        if rule != "treat":
            frozen_rows[-1].update(
                bulk_padj=b.get("bulk_padj", ""), established_state_rule=STATE_RULES[rule]
            )
        if with_untestable:
            frozen_rows[-1]["genetic_untestable_shared_posterior"] = bool_text(
                bool(g.get("genetic_untestable", False))
            )

    observed = {
        f"{rule}_rows": state_rows,
        f"{rule}_symbols": sum(
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
    if with_untestable:
        observed["genetic_untestable_named"] = sum(
            1 for value in genetics.values() if value.get("genetic_untestable", False)
        )
    names = metric_names(rule, with_untestable)
    if set(expected_counts) != set(names):
        raise ContractError(f"expected counts must give exactly {names}")
    state_label = "TREAT FDR < 0.05" if rule == "treat" else "padj < 0.05 and |logFC| > 0.5"
    audit_rows = []
    failures = []
    for metric in names:
        expected = expected_counts[metric]
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
                    "canonical_rows": "bulk rows with padj < 0.05 and |logFC| > 0.5",
                    "canonical_symbols": "unique nonempty symbols with any padj < 0.05 and |logFC| > 0.5 row",
                    "primary_susie_named": "unique named genes with SuSiE PP.H4 > 0.5",
                    "primary_susie_bulk_testable": "primary genes with a nonempty canonical bulk symbol",
                    "primary_overlap": f"primary bulk-testable genes with {state_label}"
                    if rule != "treat"
                    else "primary bulk-testable genes with TREAT FDR < 0.05",
                    "genetic_untestable_named": (
                        f"named genes without SuSiE PP.H4 > 0.5 and with a Tier-1/2 {UNTESTABLE_METHOD} pair"
                    ),
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
    if has_abf_driver:
        fields.extend(ABF_DRIVER_FIELDS)
    if rule != "treat":
        fields.extend(["bulk_padj", "established_state_rule"])
    if with_untestable:
        fields.append("genetic_untestable_shared_posterior")
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
    expected_counts = json.loads(args.expected.read_text(encoding="utf-8"))
    classes, audit = build_frozen_classes(
        sources, candidate, expected_counts, args.state_rule, args.untestable_state
    )
    write_stage_seal(
        candidate,
        "01_freeze_and_rederive",
        [manifest, classes, audit],
        [args.contract.resolve(), args.expected.resolve()],
    )
    print(f"PASS: frozen contract written beneath {candidate}")


if __name__ == "__main__":
    main()
