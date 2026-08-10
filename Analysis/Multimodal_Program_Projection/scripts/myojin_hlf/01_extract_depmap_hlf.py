#!/usr/bin/env python3
"""Extract version-matched, outcome-independent HLF covariates from DepMap."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from myojin_firewall_common import (
    CANDIDATE_ROOT,
    HLF_MODEL_ID,
    RELEASE_ID,
    atomic_write_text,
    md5_file,
    normalized_model_name,
    require_within,
    sha256_file,
    split_gene_header,
    stable_bundle_sha256,
    write_tsv,
)


PINNED_RELEASE = "DepMap Public 24Q4"
PINNED_DOI = "10.25452/figshare.plus.27993248.v1"
PINNED_FILES = {
    "Model.csv": {
        "figshare_file_id": "51065297",
        "bytes": 645696,
        "md5": "675210d17675f3517b0ce39a3c274f16",
    },
    "CRISPRGeneEffect.csv": {
        "figshare_file_id": "51064667",
        "bytes": 428678699,
        "md5": "6edf7ade09b9b34199210b559d4745d3",
    },
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": {
        "figshare_file_id": "51065489",
        "bytes": 506628654,
        "md5": "71794802b750ce77c422dad0720a40af",
    },
}


def verify_source(path: Path) -> None:
    expected = PINNED_FILES[path.name]
    if path.stat().st_size != expected["bytes"]:
        raise ValueError(
            f"Byte-size mismatch for {path.name}: {path.stat().st_size} != {expected['bytes']}"
        )
    observed_md5 = md5_file(path)
    if observed_md5 != expected["md5"]:
        raise ValueError(f"MD5 mismatch for {path.name}: {observed_md5} != {expected['md5']}")


def read_hlf_model(model_path: Path) -> dict[str, str]:
    with model_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        model_hits = []
        name_hits = []
        hlfa_hits = []
        for row in reader:
            if row.get("ModelID") == HLF_MODEL_ID:
                model_hits.append(row)
            names = [row.get("CellLineName", ""), row.get("StrippedCellLineName", "")]
            normalized = {normalized_model_name(value) for value in names if value}
            if "HLF" in normalized:
                name_hits.append(row)
            if "HLFA" in normalized:
                hlfa_hits.append(row)
    if len(model_hits) != 1:
        raise ValueError(f"Expected one {HLF_MODEL_ID} Model row; observed {len(model_hits)}")
    model = model_hits[0]
    if model.get("StrippedCellLineName") != "HLF" or model.get("CellLineName") != "HLF":
        raise ValueError(f"Model {HLF_MODEL_ID} does not map exactly to HLF")
    if model.get("OncotreeLineage") != "Liver":
        raise ValueError(f"Model {HLF_MODEL_ID} is not annotated as liver")
    if len(name_hits) != 1 or name_hits[0].get("ModelID") != HLF_MODEL_ID:
        raise ValueError("HLF name is not a unique mapping to ACH-000393")
    if any(row.get("ModelID") == HLF_MODEL_ID for row in hlfa_hits):
        raise ValueError("HLF-a was conflated with ACH-000393")
    return model


def read_wide_hlf_row(path: Path) -> tuple[list[str], list[str], str, str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if "ModelID" in header:
            model_index = header.index("ModelID")
        elif header and header[0] == "":
            model_index = 0
        else:
            model_index = 0
        default_index = next(
            (
                header.index(name)
                for name in ("IsDefaultEntryForModel", "is_default_entry")
                if name in header
            ),
            None,
        )
        hits: list[list[str]] = []
        for row in reader:
            if len(row) <= model_index:
                continue
            if row[model_index] == HLF_MODEL_ID:
                hits.append(row)
    if default_index is not None and len(hits) > 1:
        hits = [row for row in hits if row[default_index].strip().lower() in {"1", "true", "t"}]
    if len(hits) != 1:
        raise ValueError(f"Expected one default HLF row in {path.name}; observed {len(hits)}")
    row = hits[0]
    profile_id = row[header.index("ProfileID")] if "ProfileID" in header else HLF_MODEL_ID
    default_value = row[default_index] if default_index is not None else "not_applicable_single_row"
    return header, row, profile_id, default_value


def gene_value_map(header: list[str], row: list[str]) -> dict[str, tuple[str, str, str]]:
    metadata = {
        "",
        "ModelID",
        "ProfileID",
        "IsDefaultEntryForModel",
        "is_default_entry",
        "SequencingID",
    }
    result: dict[str, tuple[str, str, str]] = {}
    for column, value in zip(header, row):
        if column in metadata:
            continue
        symbol, entrez = split_gene_header(column)
        if not symbol:
            continue
        if symbol in result:
            raise ValueError(f"Duplicate normalized DepMap symbol: {symbol}")
        result[symbol] = (column, entrez, value)
    return result


def as_finite(value: str, field: str, gene: str) -> float | None:
    if value == "":
        return None
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite {field} for {gene}")
    return parsed


def merge_covariate_maps(
    effect: dict[str, tuple[str, str, str]],
    expression: dict[str, tuple[str, str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Join HLF covariates without ever resolving an Entrez conflict by symbol.

    A symbol present in both pinned DepMap matrices with different nonempty
    Entrez IDs is retained in the audit but both covariate values are blanked.
    This makes the gene untestable rather than silently assigning one gene's
    value to another identifier.
    """

    covariate_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    for symbol in sorted(set(effect) | set(expression)):
        effect_header_name, effect_entrez, effect_value = effect.get(symbol, ("", "", ""))
        expression_header_name, expression_entrez, expression_value = expression.get(
            symbol, ("", "", "")
        )
        entrez_conflict = bool(
            effect_entrez and expression_entrez and effect_entrez != expression_entrez
        )
        if entrez_conflict:
            mapping_state = "entrez_conflict_values_blanked"
            mapping_reason = (
                f"Chronos_Entrez_{effect_entrez}_differs_from_"
                f"expression_Entrez_{expression_entrez}"
            )
            chronos = None
            log2_tpm = None
            depmap_entrez = ""
        else:
            chronos = as_finite(effect_value, "Chronos", symbol)
            log2_tpm = as_finite(expression_value, "expression", symbol)
            depmap_entrez = effect_entrez or expression_entrez
            if effect_header_name and expression_header_name:
                mapping_state = (
                    "symbol_and_entrez_concordant"
                    if effect_entrez and expression_entrez
                    else "symbol_concordant_entrez_incomplete"
                )
                mapping_reason = "pass"
            elif expression_header_name:
                mapping_state = "expression_only_missing_Chronos"
                mapping_reason = "symbol_absent_from_CRISPRGeneEffect_header"
            else:
                mapping_state = "Chronos_only_missing_expression"
                mapping_reason = "symbol_absent_from_expression_header"
        tpm = None if log2_tpm is None else max(0.0, (2.0**log2_tpm) - 1.0)
        covariate_rows.append(
            {
                "gene_symbol": symbol,
                "depmap_entrez_id": depmap_entrez,
                "expression_entrez_id": expression_entrez,
                "chronos_entrez_id": effect_entrez,
                "depmap_mapping_state": mapping_state,
                "depmap_entrez_conflict": str(entrez_conflict).upper(),
                "expression_source_header": expression_header_name,
                "chronos_source_header": effect_header_name,
                "HLF_expression_log2_tpm_plus_1": (
                    "" if log2_tpm is None else f"{log2_tpm:.17g}"
                ),
                "HLF_TPM": "" if tpm is None else f"{tpm:.17g}",
                "log1p_HLF_TPM": (
                    "" if log2_tpm is None else f"{log2_tpm * math.log(2):.17g}"
                ),
                "HLF_Chronos": "" if chronos is None else f"{chronos:.17g}",
                "expression_mapped": str(log2_tpm is not None).upper(),
                "chronos_mapped": str(chronos is not None).upper(),
                "depmap_release": PINNED_RELEASE,
                "model_id": HLF_MODEL_ID,
            }
        )
        audit_rows.append(
            {
                "gene_symbol": symbol,
                "expression_source_header": expression_header_name,
                "expression_entrez_id": expression_entrez,
                "chronos_source_header": effect_header_name,
                "chronos_entrez_id": effect_entrez,
                "depmap_mapping_state": mapping_state,
                "depmap_entrez_conflict": str(entrez_conflict).upper(),
                "covariate_values_blanked": str(entrez_conflict).upper(),
                "mapping_reason": mapping_reason,
                "eligible_for_primary_mapping": str(
                    not entrez_conflict and bool(expression_header_name and effect_header_name)
                ).upper(),
                "release_id": RELEASE_ID,
            }
        )
    return covariate_rows, audit_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=CANDIDATE_ROOT)
    parser.add_argument("--depmap-release", default=PINNED_RELEASE)
    args = parser.parse_args()

    outdir = require_within(args.outdir, CANDIDATE_ROOT)
    source_dir = require_within(args.source_dir, outdir)
    if args.depmap_release != PINNED_RELEASE:
        raise ValueError(
            f"This frozen adapter requires {PINNED_RELEASE}; observed {args.depmap_release}"
        )
    paths = {name: source_dir / name for name in PINNED_FILES}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
        verify_source(path)

    model = read_hlf_model(paths["Model.csv"])
    effect_header, effect_row, effect_profile, effect_default = read_wide_hlf_row(
        paths["CRISPRGeneEffect.csv"]
    )
    expression_header, expression_row, expression_profile, expression_default = read_wide_hlf_row(
        paths["OmicsExpressionProteinCodingGenesTPMLogp1.csv"]
    )
    effect = gene_value_map(effect_header, effect_row)
    expression = gene_value_map(expression_header, expression_row)
    covariate_rows, mapping_audit_rows = merge_covariate_maps(effect, expression)

    covariates = outdir / "detectability_covariates.tsv"
    write_tsv(
        covariates,
        covariate_rows,
        [
            "gene_symbol",
            "depmap_entrez_id",
            "expression_entrez_id",
            "chronos_entrez_id",
            "depmap_mapping_state",
            "depmap_entrez_conflict",
            "expression_source_header",
            "chronos_source_header",
            "HLF_expression_log2_tpm_plus_1",
            "HLF_TPM",
            "log1p_HLF_TPM",
            "HLF_Chronos",
            "expression_mapped",
            "chronos_mapped",
            "depmap_release",
            "model_id",
        ],
    )

    mapping_audit = outdir / "depmap_gene_mapping_audit.tsv"
    write_tsv(
        mapping_audit,
        mapping_audit_rows,
        [
            "gene_symbol",
            "expression_source_header",
            "expression_entrez_id",
            "chronos_source_header",
            "chronos_entrez_id",
            "depmap_mapping_state",
            "depmap_entrez_conflict",
            "covariate_values_blanked",
            "mapping_reason",
            "eligible_for_primary_mapping",
            "release_id",
        ],
    )

    identity = outdir / "hlf_identity_audit.tsv"
    write_tsv(
        identity,
        [
            {
                "model_id": model.get("ModelID", ""),
                "cell_line_name": model.get("CellLineName", ""),
                "stripped_cell_line_name": model.get("StrippedCellLineName", ""),
                "depmap_model_type": model.get("DepmapModelType", ""),
                "oncotree_lineage": model.get("OncotreeLineage", ""),
                "oncotree_primary_disease": model.get("OncotreePrimaryDisease", ""),
                "ccle_name": model.get("CCLEName", ""),
                "rrid": model.get("RRID", ""),
                "expression_profile_id": expression_profile,
                "expression_default_entry": expression_default,
                "chronos_profile_id": effect_profile,
                "chronos_default_entry": effect_default,
                "hlf_unique": "TRUE",
                "hlfa_conflated": "FALSE",
                "actual_vehicle_fitness_observed": "FALSE",
                "chronos_role": "generic_HLF_fitness_covariate_not_Myojin_vehicle_effect",
            }
        ],
        [
            "model_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "depmap_model_type",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "ccle_name",
            "rrid",
            "expression_profile_id",
            "expression_default_entry",
            "chronos_profile_id",
            "chronos_default_entry",
            "hlf_unique",
            "hlfa_conflated",
            "actual_vehicle_fitness_observed",
            "chronos_role",
        ],
    )

    depmap_manifest = outdir / "depmap_source_manifest.tsv"
    manifest_rows = []
    for name, path in paths.items():
        metadata = PINNED_FILES[name]
        manifest_rows.append(
            {
                "depmap_release": PINNED_RELEASE,
                "doi": PINNED_DOI,
                "figshare_article_id": "27993248",
                "figshare_file_id": metadata["figshare_file_id"],
                "filename": name,
                "bytes": path.stat().st_size,
                "published_md5": metadata["md5"],
                "observed_md5": md5_file(path),
                "observed_sha256": sha256_file(path),
                "license": "CC_BY_4.0",
                "release_id": RELEASE_ID,
            }
        )
    write_tsv(
        depmap_manifest,
        manifest_rows,
        [
            "depmap_release",
            "doi",
            "figshare_article_id",
            "figshare_file_id",
            "filename",
            "bytes",
            "published_md5",
            "observed_md5",
            "observed_sha256",
            "license",
            "release_id",
        ],
    )
    conflicts = sum(
        row["depmap_entrez_conflict"] == "TRUE" for row in mapping_audit_rows
    )
    bundle_hash, _ = stable_bundle_sha256(
        [covariates, mapping_audit, identity, depmap_manifest], outdir
    )
    atomic_write_text(
        outdir / "DEPMAP_READY",
        "release_id\tstatus\tdepmap_release\tbundle_sha256\tentrez_conflicts_excluded\toutcomes_read\n"
        f"{RELEASE_ID}\tdepmap_hlf_covariates_ready\t{PINNED_RELEASE}\t{bundle_hash}\t{conflicts}\tFALSE\n",
    )
    print(
        f"DEPMAP_READY genes={len(covariate_rows)} model={HLF_MODEL_ID} "
        f"release={PINNED_RELEASE!r} entrez_conflicts_excluded={conflicts} "
        f"outcomes_read=FALSE bundle_sha256={bundle_hash}"
    )


if __name__ == "__main__":
    main()
