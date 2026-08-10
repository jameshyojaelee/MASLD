#!/usr/bin/env python3
"""Fail closed unless the active rendered figure-label audit matches its PDFs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import re
import subprocess
from collections import Counter
from pathlib import Path, PurePosixPath


SCOPE_FIELDS = ("panel_id", "pdf_path", "label_mode", "classification_reason")
ROSTER_FIELDS = ("panel_id", "producer_script", "pdf_path", "scope_reason")
INVENTORY_FIELDS = SCOPE_FIELDS + (
    "pdf_sha256",
    "pdftotext_version",
    "extracted_text_sha256",
)
TOKEN_FIELDS = (
    "panel_id",
    "pdf_path",
    "token",
    "classification",
    "reason",
    "proof_source_path",
    "proof_source_line",
)
PROTECTED_FIELDS = ("symbol", "source_path", "source_line", "context", "baseline_status")
VALID_MODES = {"dynamic_rendered", "literal_static", "no_gene_labels"}
EXPECTED_MODE_COUNTS = Counter(
    {"dynamic_rendered": 24, "literal_static": 9, "no_gene_labels": 44}
)
TOKEN_RE = re.compile(
    r"(?<![A-Za-z0-9_.-])([A-Za-z0-9][A-Za-z0-9_.-]*)(?![A-Za-z0-9_.-])"
)
DYNAMIC_GENE_REASON = (
    "Data-selected or data-filtered rendered gene label; protected from the frozen PDF."
)
STATIC_GENE_REASON = (
    "Static or literal rendered gene label; source-proven by protected_figure_genes.tsv."
)


def _pip() -> str:
    return "GENCODE collision: PIP denotes posterior inclusion probability, not the PIP gene."


def _stage(token: str) -> str:
    return f"GENCODE collision: {token} denotes a fibrosis stage in this panel, not the {token} gene."


# Exact panel-token pairs are deliberately allowlisted. This prevents a real
# rendered gene from being reclassified as non-gene merely to remove it from
# the protected set.
EXPECTED_NON_GENE_COLLISIONS = {
    ("fig2b_pip_vs_susie_coloc", "PIP"): _pip(),
    ("fig2f_rora_locus", "PIP"): _pip(),
    ("fig2g_fabp1_locus", "PIP"): _pip(),
    ("fig2h_locus_legend", "CS"): (
        "GENCODE collision: CS denotes credible set in the locus legend, not the CS gene."
    ),
    ("figs2a_gwas_perstudy", "MVP"): (
        "GENCODE collision: MVP denotes the Million Veteran Program, not the MVP gene."
    ),
    ("figs2b_finemap_cascade", "PIP"): _pip(),
    ("figs2c_finemap_headtohead", "PIP"): _pip(),
    ("figs2d_efhd1_locus", "CS"): (
        "GENCODE collision: CS denotes credible set in fine-mapping text, not the CS gene."
    ),
    ("figs2d_efhd1_locus", "PIP"): _pip(),
    ("figs2e_acads_locus", "CS"): (
        "GENCODE collision: CS denotes credible set in fine-mapping text, not the CS gene."
    ),
    ("figs2e_acads_locus", "PIP"): _pip(),
    ("figs2f_coloc_concordance", "PIP"): _pip(),
    ("figs2j_pip_concentration", "PIP"): _pip(),
    ("figs2q_scchromatin", "GRN"): (
        "GENCODE collision: GRN denotes gene-regulatory network in the annotation text, not the GRN gene."
    ),
    ("fig3b_nas_fib_grid", "F2"): _stage("F2"),
    ("fig3b_nas_fib_grid", "F3"): _stage("F3"),
    ("fig3e_progression", "F2"): _stage("F2"),
    ("fig3e_progression", "F3"): _stage("F3"),
    ("fig3i_tf_convergence", "F2"): _stage("F2"),
    ("fig3i_tf_convergence", "F3"): _stage("F3"),
    ("fig3i_tf_convergence", "TF"): (
        "GENCODE collision: TF denotes transcription factor in an axis label, not the TF gene."
    ),
    ("figs3g_fib_degs", "F2"): _stage("F2"),
    ("figs3g_fib_degs", "F3"): _stage("F3"),
    ("figs3h_fib_upset", "F2"): _stage("F2"),
    ("figs3h_fib_upset", "F3"): _stage("F3"),
    ("figs3i_deg_heterogeneity", "C2"): (
        "GENCODE collision: C2 denotes the C2 statistical model in subtitle text, not complement C2."
    ),
    ("figs3k_power_curves", "TPR"): (
        "GENCODE collision: TPR denotes true-positive rate in an axis label, not the TPR gene."
    ),
    ("figs3k_simulation_power", "TPR"): (
        "GENCODE collision: TPR denotes true-positive rate in an axis label, not the TPR gene."
    ),
    ("figs3o_fib_direction", "F2"): _stage("F2"),
    ("figs3o_fib_direction", "F3"): _stage("F3"),
    ("figs3p_hallmark", "F2"): _stage("F2"),
    ("figs3p_hallmark", "F3"): _stage("F3"),
    ("figs3p_hallmark", "NES"): (
        "GENCODE collision: NES denotes normalized enrichment score, not the NES gene."
    ),
    ("figs3s_carrier_routing", "F2"): _stage("F2"),
    ("figs3s_carrier_routing", "F3"): _stage("F3"),
    ("figs3x_geneset_correlation", "F3"): (
        "GENCODE collision: F3 is a rendered stage/text fragment in this panel, not the F3 gene."
    ),
    ("figs3x_geneset_correlation", "FH"): (
        "GENCODE collision: FH is a rotated pathway-label text fragment, not the FH gene."
    ),
    ("fig4c_mrna_protein", "NES"): (
        "GENCODE collision: NES denotes normalized enrichment score in an axis label, not the NES gene."
    ),
    ("fig5a_therapeutic_axes", "PIP"): _pip(),
    ("fig5e_glp1ra_axis", "C2"): (
        "GENCODE collision: C2 is a rendered cell-class label fragment, not complement C2."
    ),
}


class InventoryError(RuntimeError):
    """Raised when a rendered-label binding or classification gate fails."""


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_exact_tsv(path: Path, expected_fields: tuple[str, ...]) -> list[dict[str, str]]:
    if not path.is_file():
        raise InventoryError(f"Missing figure-label sidecar: {path}")
    payload = path.read_bytes()
    if b"\r" in payload:
        raise InventoryError(f"CR/CRLF is forbidden in machine-readable sidecar: {path}")
    try:
        text = payload.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise InventoryError(f"Sidecar is not valid UTF-8: {path}") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
    if tuple(reader.fieldnames or ()) != expected_fields:
        raise InventoryError(
            f"Unexpected columns in {path}: {reader.fieldnames}; expected {expected_fields}"
        )
    rows = list(reader)
    if not rows:
        raise InventoryError(f"Empty figure-label sidecar: {path}")
    for line_number, row in enumerate(rows, start=2):
        for value in row.values():
            if value is None or not value.strip() or "\n" in value or "\r" in value or "\t" in value:
                raise InventoryError(f"Blank or embedded-control field in {path}:{line_number}")
    return rows


def checked_relative(value: str, *, prefix: str, suffix: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != value:
        raise InventoryError(f"Unsafe or noncanonical relative path: {value}")
    if not value.startswith(prefix) or not value.endswith(suffix):
        raise InventoryError(f"Path outside figure-label contract: {value}")
    return path


def resolve_inside(project_root: Path, relative: PurePosixPath) -> Path:
    candidate = project_root.joinpath(*relative.parts)
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(project_root)
    except (FileNotFoundError, ValueError) as exc:
        raise InventoryError(f"Missing or escaping figure-label path: {relative}") from exc
    if not resolved.is_file():
        raise InventoryError(f"Figure-label path is not a file: {relative}")
    return resolved


def load_scope(project_root: Path, scope_path: Path) -> list[dict[str, str]]:
    rows = read_exact_tsv(scope_path, SCOPE_FIELDS)
    panel_ids: set[str] = set()
    pdf_paths: set[str] = set()
    modes: Counter[str] = Counter()
    for row in rows:
        panel_id = row["panel_id"]
        if not re.fullmatch(r"[a-z0-9_]+", panel_id) or panel_id in panel_ids:
            raise InventoryError(f"Invalid or duplicate panel_id: {panel_id}")
        if row["pdf_path"] in pdf_paths:
            raise InventoryError(f"Duplicate active PDF: {row['pdf_path']}")
        if row["label_mode"] not in VALID_MODES:
            raise InventoryError(f"Unclassified active panel {panel_id}: {row['label_mode']}")
        panel_ids.add(panel_id)
        pdf_paths.add(row["pdf_path"])
        modes[row["label_mode"]] += 1
        pdf = checked_relative(row["pdf_path"], prefix="figures/", suffix=".pdf")
        resolve_inside(project_root, pdf)
    if modes != EXPECTED_MODE_COUNTS:
        raise InventoryError(
            f"Active-panel classification changed: actual={dict(modes)} expected={dict(EXPECTED_MODE_COUNTS)}"
        )
    if "figs3_module_heatmap_full" not in panel_ids:
        raise InventoryError("Active full 54-module supplement is missing from the panel scope")
    collision_panels = {panel for panel, _ in EXPECTED_NON_GENE_COLLISIONS}
    if not collision_panels.issubset(panel_ids):
        raise InventoryError(
            "Collision allowlist refers to panels absent from active scope: "
            + ",".join(sorted(collision_panels - panel_ids))
        )
    return rows


def load_dynamic_roster(
    project_root: Path,
    roster_path: Path,
    scope_by_panel: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows = read_exact_tsv(roster_path, ROSTER_FIELDS)
    roster_ids: set[str] = set()
    for row in rows:
        panel_id = row["panel_id"]
        if panel_id in roster_ids or panel_id not in scope_by_panel:
            raise InventoryError(f"Duplicate or unscoped dynamic panel: {panel_id}")
        scope = scope_by_panel[panel_id]
        if scope["label_mode"] != "dynamic_rendered" or scope["pdf_path"] != row["pdf_path"]:
            raise InventoryError(f"Dynamic roster/scope disagreement for {panel_id}")
        producer = checked_relative(row["producer_script"], prefix="scripts/figures/", suffix=".R")
        resolve_inside(project_root, producer)
        roster_ids.add(panel_id)
    expected_ids = {
        panel_id
        for panel_id, row in scope_by_panel.items()
        if row["label_mode"] == "dynamic_rendered"
    }
    if roster_ids != expected_ids:
        raise InventoryError(
            "Dynamic rendered roster does not cover active scope exactly: "
            f"missing={sorted(expected_ids - roster_ids)} extra={sorted(roster_ids - expected_ids)}"
        )
    return rows


def load_protected_static(path: Path) -> dict[str, dict[str, str]]:
    rows = read_exact_tsv(path, PROTECTED_FIELDS)
    by_symbol: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["symbol"] in by_symbol:
            raise InventoryError(f"Duplicate symbol in static protected inventory: {row['symbol']}")
        by_symbol[row["symbol"]] = row
    return by_symbol


def pdftotext_version() -> str:
    try:
        result = subprocess.run(
            ["pdftotext", "-v"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise InventoryError("pdftotext is unavailable or failed its version check") from exc
    lines = result.stdout.decode("utf-8", errors="strict").splitlines()
    if not lines or not lines[0].startswith("pdftotext version "):
        raise InventoryError(f"Unrecognized pdftotext version output: {lines!r}")
    return lines[0]


def extract_pdf_text(pdf: Path) -> bytes:
    try:
        result = subprocess.run(
            ["pdftotext", "-layout", str(pdf), "-"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        detail = ""
        if isinstance(exc, subprocess.CalledProcessError) and exc.stderr:
            detail = ": " + exc.stderr.decode("utf-8", errors="replace").strip()
        raise InventoryError(f"pdftotext failed for {pdf}{detail}") from exc
    result.stdout.decode("utf-8", errors="strict")
    return result.stdout


def official_symbols(project_root: Path) -> set[str]:
    metadata = project_root / "data/gencode_v49_gene_metadata.tsv.gz"
    if not metadata.is_file():
        raise InventoryError(f"Missing GENCODE v49 metadata: {metadata}")
    with gzip.open(metadata, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "gene_name" not in reader.fieldnames:
            raise InventoryError(f"GENCODE metadata lacks gene_name: {metadata}")
        symbols = {row["gene_name"].strip() for row in reader if row["gene_name"].strip()}
    if not symbols:
        raise InventoryError(f"No official symbols loaded from {metadata}")
    return symbols


def extract_official_tokens(text: bytes, official: set[str]) -> set[str]:
    rendered = text.decode("utf-8", errors="strict")
    return {token for token in TOKEN_RE.findall(rendered) if token in official}


def inspect_scope(
    project_root: Path,
    scope_rows: list[dict[str, str]],
) -> tuple[str, list[dict[str, str]], dict[str, set[str]]]:
    version = pdftotext_version()
    official = official_symbols(project_root)
    inventory_rows: list[dict[str, str]] = []
    tokens_by_panel: dict[str, set[str]] = {}
    for row in scope_rows:
        pdf = resolve_inside(
            project_root,
            checked_relative(row["pdf_path"], prefix="figures/", suffix=".pdf"),
        )
        text = extract_pdf_text(pdf)
        inventory_rows.append(
            {
                **row,
                "pdf_sha256": sha256_file(pdf),
                "pdftotext_version": version,
                "extracted_text_sha256": sha256_bytes(text),
            }
        )
        tokens_by_panel[row["panel_id"]] = extract_official_tokens(text, official)
    return version, inventory_rows, tokens_by_panel


def expected_token_row(
    panel: dict[str, str],
    token: str,
    protected_static: dict[str, dict[str, str]],
) -> dict[str, str]:
    panel_id = panel["panel_id"]
    collision_reason = EXPECTED_NON_GENE_COLLISIONS.get((panel_id, token))
    if collision_reason:
        return {
            "panel_id": panel_id,
            "pdf_path": panel["pdf_path"],
            "token": token,
            "classification": "non_gene_collision",
            "reason": collision_reason,
            "proof_source_path": panel["pdf_path"],
            "proof_source_line": "rendered_context",
        }
    mode = panel["label_mode"]
    if mode == "no_gene_labels":
        raise InventoryError(
            f"Panel classified no_gene_labels contains unapproved official-symbol token: {(panel_id, token)}"
        )
    if mode == "literal_static":
        proof = protected_static.get(token)
        if proof is None:
            raise InventoryError(
                f"Static rendered gene lacks source-level protected proof: {(panel_id, token)}"
            )
        return {
            "panel_id": panel_id,
            "pdf_path": panel["pdf_path"],
            "token": token,
            "classification": "gene_label",
            "reason": STATIC_GENE_REASON,
            "proof_source_path": proof["source_path"],
            "proof_source_line": proof["source_line"],
        }
    return {
        "panel_id": panel_id,
        "pdf_path": panel["pdf_path"],
        "token": token,
        "classification": "gene_label",
        "reason": DYNAMIC_GENE_REASON,
        "proof_source_path": panel["pdf_path"],
        "proof_source_line": "rendered_pdf",
    }


def validate_inventory(
    project_root: Path,
    scope_path: Path,
    roster_path: Path,
    inventory_path: Path,
    token_path: Path,
    protected_path: Path,
) -> dict[str, int | str]:
    project_root = project_root.resolve(strict=True)
    scope_rows = load_scope(project_root, scope_path)
    scope_by_panel = {row["panel_id"]: row for row in scope_rows}
    load_dynamic_roster(project_root, roster_path, scope_by_panel)
    protected_static = load_protected_static(protected_path)
    inventory_rows = read_exact_tsv(inventory_path, INVENTORY_FIELDS)
    token_rows = read_exact_tsv(token_path, TOKEN_FIELDS)

    if len(inventory_rows) != len(scope_rows):
        raise InventoryError("Active figure inventory and scope row counts differ")
    inventory_by_panel: dict[str, dict[str, str]] = {}
    for row in inventory_rows:
        panel_id = row["panel_id"]
        if panel_id in inventory_by_panel or panel_id not in scope_by_panel:
            raise InventoryError(f"Duplicate or unscoped inventory panel: {panel_id}")
        if any(row[field] != scope_by_panel[panel_id][field] for field in SCOPE_FIELDS):
            raise InventoryError(f"Inventory classification/path differs from scope for {panel_id}")
        for field in ("pdf_sha256", "extracted_text_sha256"):
            if not re.fullmatch(r"[0-9a-f]{64}", row[field]):
                raise InventoryError(f"Invalid {field} for {panel_id}")
        inventory_by_panel[panel_id] = row
    if set(inventory_by_panel) != set(scope_by_panel):
        raise InventoryError("Active figure inventory does not cover scope exactly")

    actual_version, actual_inventory, actual_tokens = inspect_scope(project_root, scope_rows)
    if {row["pdftotext_version"] for row in inventory_rows} != {actual_version}:
        raise InventoryError(
            f"pdftotext version differs from frozen inventory: actual={actual_version!r}"
        )
    for actual in actual_inventory:
        expected = inventory_by_panel[actual["panel_id"]]
        if actual != expected:
            differing = [field for field in INVENTORY_FIELDS if actual[field] != expected[field]]
            raise InventoryError(
                f"Active rendered PDF/text binding changed for {actual['panel_id']}: {','.join(differing)}"
            )

    expected_rows: dict[tuple[str, str], dict[str, str]] = {}
    for panel_id, panel in scope_by_panel.items():
        for token in actual_tokens[panel_id]:
            row = expected_token_row(panel, token, protected_static)
            expected_rows[(panel_id, token)] = row
    collision_keys = {
        key for key, row in expected_rows.items() if row["classification"] == "non_gene_collision"
    }
    if collision_keys != set(EXPECTED_NON_GENE_COLLISIONS):
        raise InventoryError(
            "Exact non-gene collision allowlist differs from rendered tokens: "
            f"missing={sorted(set(EXPECTED_NON_GENE_COLLISIONS) - collision_keys)} "
            f"extra={sorted(collision_keys - set(EXPECTED_NON_GENE_COLLISIONS))}"
        )

    observed_rows: dict[tuple[str, str], dict[str, str]] = {}
    for row in token_rows:
        key = (row["panel_id"], row["token"])
        if key in observed_rows:
            raise InventoryError(f"Duplicate official-symbol token classification: {key}")
        if key not in expected_rows:
            raise InventoryError(f"Extra or unrendered official-symbol token classification: {key}")
        if row != expected_rows[key]:
            differing = [field for field in TOKEN_FIELDS if row[field] != expected_rows[key][field]]
            raise InventoryError(f"Unapproved token classification for {key}: {','.join(differing)}")
        observed_rows[key] = row
    missing_rows = set(expected_rows) - set(observed_rows)
    if missing_rows:
        raise InventoryError(
            "Rendered token classification mismatch; missing=" + repr(sorted(missing_rows))
        )

    gene_rows = [row for row in token_rows if row["classification"] == "gene_label"]
    dynamic_gene_rows = [
        row
        for row in gene_rows
        if scope_by_panel[row["panel_id"]]["label_mode"] == "dynamic_rendered"
    ]
    unique_dynamic = {row["token"] for row in dynamic_gene_rows}
    return {
        "panels": len(scope_rows),
        "dynamic_panels": EXPECTED_MODE_COUNTS["dynamic_rendered"],
        "static_panels": EXPECTED_MODE_COUNTS["literal_static"],
        "no_gene_panels": EXPECTED_MODE_COUNTS["no_gene_labels"],
        "classified_rows": len(token_rows),
        "gene_label_rows": len(gene_rows),
        "dynamic_gene_label_rows": len(dynamic_gene_rows),
        "unique_dynamic_gene_labels": len(unique_dynamic),
        "non_gene_collision_rows": len(collision_keys),
        "pdftotext_version": actual_version,
    }


def default_paths(project_root: Path) -> tuple[Path, Path, Path, Path, Path]:
    base = project_root / "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
    return (
        base / "active_figure_panel_scope.tsv",
        base / "rendered_figure_label_roster.tsv",
        base / "active_figure_panel_inventory.tsv",
        base / "rendered_figure_label_tokens.tsv",
        base / "protected_figure_genes.tsv",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--scope", type=Path)
    parser.add_argument("--roster", type=Path)
    parser.add_argument("--inventory", type=Path)
    parser.add_argument("--tokens", type=Path)
    parser.add_argument("--protected", type=Path)
    args = parser.parse_args()
    project_root = args.project_root.resolve(strict=True)
    defaults = default_paths(project_root)
    try:
        result = validate_inventory(
            project_root,
            args.scope or defaults[0],
            args.roster or defaults[1],
            args.inventory or defaults[2],
            args.tokens or defaults[3],
            args.protected or defaults[4],
        )
    except (InventoryError, UnicodeDecodeError) as exc:
        raise SystemExit(f"FAIL rendered figure-label inventory: {exc}") from exc
    print(
        "PASS rendered figure-label inventory: "
        f"{result['panels']} PDFs "
        f"({result['dynamic_panels']} dynamic, {result['static_panels']} static, "
        f"{result['no_gene_panels']} no-gene); "
        f"{result['classified_rows']} official-symbol rows; "
        f"{result['unique_dynamic_gene_labels']} unique dynamic gene labels; "
        f"{result['non_gene_collision_rows']} documented collision rows; "
        f"{result['pdftotext_version']}"
    )


if __name__ == "__main__":
    main()
