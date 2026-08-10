#!/usr/bin/env python3
"""Render deterministic candidate-only individual PDF panels from frozen source tables."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from rel01_snapshot_candidate import (
    materialize_immutable,
    materialize_immutable_json,
    materialize_immutable_tsv,
)
from rel05_validate_candidate import validate_snapshot_base
from release_common import (
    CANDIDATE_ID,
    ReleaseContractError,
    read_tsv_exact,
    resolve_project_path,
    sha256_file,
)
from release_products import (
    EXPECTED_SUPPLEMENTARY_PANELS,
    FIGURE_SOURCE_MANIFEST_FIELDS,
    PANEL_BLUEPRINT_FIELDS,
    SOURCE_ROW_FIELDS,
    candidate_figure_root,
    ensure_candidate_output_root,
    load_frozen_input_index,
)


def pdf_escape(value: str) -> str:
    safe = value.encode("ascii", "replace").decode("ascii")
    return safe.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def build_pdf(panel_title: str, rows: list[dict[str, str]]) -> bytes:
    marks = [row for row in rows if row["plot_role"] == "mark"]
    if not marks:
        raise ReleaseContractError(f"panel has no plot marks: {panel_title}")
    annotations: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["plot_role"] == "annotation":
            annotations[row["group_id"]].append(row)
    non_effect_states = {"untestable", "skipped", "not_applicable"}
    quantitative_marks = [
        row for row in marks if row.get("evidence_status") not in non_effect_states
    ]
    values = [float(row["value"]) for row in quantitative_marks] or [0.0]
    scale = max(max(abs(value) for value in values), 1e-12)
    width = 600.0
    height = 400.0
    left = 45.0
    right = 575.0
    baseline = 185.0
    usable = 120.0
    slot = (right - left) / max(len(marks), 1)
    bar_width = min(28.0, slot * 0.55)
    commands = [
        "0 0 0 rg",
        f"BT /F1 6 Tf 35 375 Td ({pdf_escape(panel_title)}) Tj ET",
        "0.4 w 0 0 0 RG",
        f"{left:.2f} {baseline:.2f} m {right:.2f} {baseline:.2f} l S",
    ]
    for index, row in enumerate(marks):
        center = left + slot * (index + 0.5)
        evidence_status = row.get("evidence_status")
        if evidence_status in non_effect_states:
            box_bottom = baseline - 5.0
            glyph_commands = [
                "0.620 0.620 0.620 RG",
                f"{center - bar_width / 2:.2f} {box_bottom:.2f} "
                f"{bar_width:.2f} {bar_width:.2f} re S",
            ]
            if evidence_status in {"untestable", "not_applicable"}:
                glyph_commands.append(
                    f"{center - bar_width / 2:.2f} {box_bottom:.2f} m "
                    f"{center + bar_width / 2:.2f} {box_bottom + bar_width:.2f} l S",
                )
            if evidence_status == "untestable":
                glyph_commands.append(
                    f"{center - bar_width / 2:.2f} {box_bottom + bar_width:.2f} m "
                    f"{center + bar_width / 2:.2f} {box_bottom:.2f} l S"
                )
            glyph_commands.extend(
                [
                    "0 0 0 RG 0 0 0 rg",
                    "BT /F1 6 Tf "
                    f"{center - bar_width / 2:.2f} {baseline - 18:.2f} Td "
                    f"({pdf_escape(row['label'][:18])}) Tj ET",
                    "BT /F1 6 Tf "
                    f"{center - bar_width / 2:.2f} {box_bottom + bar_width + 5:.2f} Td "
                    f"({pdf_escape(row['display_value'])}) Tj ET",
                ]
            )
            commands.extend(glyph_commands)
            continue
        value = float(row["value"])
        extent = usable * value / scale
        bottom = baseline if extent >= 0 else baseline + extent
        bar_height = abs(extent)
        if row["is_control"] == "true":
            color = "0.620 0.620 0.620"
        elif value >= 0:
            color = "0.298 0.447 0.690"
        else:
            color = "0.769 0.306 0.322"
        commands.extend(
            [
                f"{color} rg",
                f"{center - bar_width / 2:.2f} {bottom:.2f} "
                f"{bar_width:.2f} {max(bar_height, 0.8):.2f} re f",
                "0 0 0 rg",
                "BT /F1 6 Tf "
                f"{center - bar_width / 2:.2f} {baseline - 18:.2f} Td "
                f"({pdf_escape(row['label'][:18])}) Tj ET",
                "BT /F1 6 Tf "
                f"{center - bar_width / 2:.2f} "
                f"{(bottom + bar_height + 5 if extent >= 0 else bottom - 9):.2f} Td "
                f"({pdf_escape(row['display_value'])}) Tj ET",
            ]
        )
        note_rows = sorted(
            annotations.get(row["group_id"], []), key=lambda item: item["number_role"]
        )
        for note_index, note_row in enumerate(note_rows):
            note = f"{note_row['number_role']}={note_row['display_value']}"
            commands.append(
                "BT /F1 6 Tf "
                f"{center - bar_width / 2:.2f} "
                f"{baseline - 29 - note_index * 8:.2f} Td "
                f"({pdf_escape(note)}) Tj ET"
            )
    mark_group_ids = {row["group_id"] for row in marks}
    orphan_annotations = sorted(
        [
            row
            for row in rows
            if row["plot_role"] == "annotation"
            and row["group_id"] not in mark_group_ids
        ],
        key=lambda item: (int(item.get("plot_order", "0")), item["number_role"]),
    )
    for index, row in enumerate(orphan_annotations):
        text_value = f"{row['number_role']}={row['display_value']} ({row['label']})"
        commands.append(
            "BT /F1 6 Tf "
            f"330 {355 - index * 10:.2f} Td "
            f"({pdf_escape(text_value)}) Tj ET"
        )
    commands.extend(
        [
            "BT /F1 6 Tf 35 20 Td "
            f"({pdf_escape('All prespecified rows retained; status is not a plotting filter.')}) Tj ET",
            "BT /F1 6 Tf 35 10 Td "
            f"({pdf_escape('Neutral glyphs: X=untestable, /=not applicable, hollow=skipped; none are zero effects.')}) Tj ET",
        ]
    )
    content = ("\n".join(commands) + "\n").encode("ascii")
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width:.0f} {height:.0f}] "
            "/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>"
        ).encode("ascii"),
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        b"<< /Length "
        + str(len(content)).encode("ascii")
        + b" >>\nstream\n"
        + content
        + b"endstream",
    ]
    payload = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0]
    for number, obj in enumerate(objects, start=1):
        offsets.append(len(payload))
        payload.extend(f"{number} 0 obj\n".encode("ascii"))
        payload.extend(obj)
        payload.extend(b"\nendobj\n")
    xref_offset = len(payload)
    payload.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    payload.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        payload.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
    payload.extend(
        (
            f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
            f"startxref\n{xref_offset}\n%%EOF\n"
        ).encode("ascii")
    )
    return bytes(payload)


def validate_rel02_state(
    candidate: Path, spec_hash: str, fixture_mode: bool
) -> dict[str, object]:
    state_path = candidate / "manifests/rel02_state.json"
    from release_products import read_json_object

    state = read_json_object(state_path, "REL-02 state")
    expected = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "state": "REL02_TABLES_FROZEN",
        "input_snapshot_spec_sha256": spec_hash,
        "journal_branch": "Cell Genomics Resource",
        "figure_count": 5,
        "myojin_role": "supplement_only",
        "n_supplementary_panels": sum(
            len(panels) for panels in EXPECTED_SUPPLEMENTARY_PANELS.values()
        ),
        "fixture_mode": fixture_mode,
        "canonical_promotion_status": "not_promoted",
    }
    for key, value in expected.items():
        if state.get(key) != value:
            raise ReleaseContractError(f"REL-02 state drift for {key}")
    hash_paths = {
        "numbers_ledger_sha256": candidate / "tables/numbers_ledger.tsv",
        "claim_ledger_sha256": candidate / "claims/claim_ledger.tsv",
        "source_dependency_ledger_sha256": candidate
        / "claims/source_dependency_ledger.tsv",
        "acceptance_gates_sha256": candidate / "manifests/acceptance_gates.tsv",
        "figure_panel_blueprint_sha256": candidate
        / "manifests/figure_panel_blueprint.tsv",
        "analysis_release_manifest_sha256": candidate
        / "manifests/analysis_release_manifest.tsv",
    }
    for field, path in hash_paths.items():
        if state.get(field) != sha256_file(path):
            raise ReleaseContractError(f"REL-02 product hash drift: {field}")
    return state


def render_candidate_panels(
    project_root: Path, fixture_mode: bool = False
) -> dict[str, object]:
    base = validate_snapshot_base(project_root, fixture_mode=fixture_mode)
    project = base["project"]
    candidate = base["candidate"]
    if not isinstance(project, Path) or not isinstance(candidate, Path):
        raise ReleaseContractError("internal REL-03 path type drift")
    load_frozen_input_index(project, candidate, str(base["spec_hash"]), fixture_mode)
    rel02_state = validate_rel02_state(candidate, str(base["spec_hash"]), fixture_mode)
    expected_figure_root = candidate_figure_root(project)
    blueprint_rows = read_tsv_exact(
        candidate / "manifests/figure_panel_blueprint.tsv", PANEL_BLUEPRINT_FIELDS
    )
    if not blueprint_rows:
        raise ReleaseContractError("figure panel blueprint is empty")
    expected_supplementary = {
        (figure_id, panel[0])
        for figure_id, panels in EXPECTED_SUPPLEMENTARY_PANELS.items()
        for panel in panels
    }
    observed_supplementary = {
        (row["figure_id"], row["panel_id"])
        for row in blueprint_rows
        if row["figure_id"].startswith("FigureS")
    }
    if observed_supplementary != expected_supplementary:
        raise ReleaseContractError(
            "render blueprint does not contain exactly the locked supplementary panels: "
            f"observed={sorted(observed_supplementary)}"
        )
    identities = set()
    prepared_panels = []
    for row in blueprint_rows:
        identity = (row["figure_id"], row["panel_id"])
        if identity in identities:
            raise ReleaseContractError(f"duplicate figure panel identity: {identity}")
        identities.add(identity)
        if row["normal_font_pt"] != "6" or row["format"] != "PDF":
            raise ReleaseContractError(f"panel typography/format drift: {identity}")
        source_relative = row["source_table_path"]
        source = candidate / source_relative
        if (
            not source_relative.startswith("figure_sources/")
            or not source.is_file()
            or source.is_symlink()
            or sha256_file(source) != row["source_table_sha256"]
        ):
            raise ReleaseContractError(f"candidate source table drift: {identity}")
        source_rows = read_tsv_exact(source, SOURCE_ROW_FIELDS)
        if not source_rows or any(
            item["figure_id"] != row["figure_id"] or item["panel_id"] != row["panel_id"]
            for item in source_rows
        ):
            raise ReleaseContractError(f"panel/source identity mismatch: {identity}")
        pdf_path = resolve_project_path(
            project, row["pdf_path"], f"{identity} PDF path"
        )
        try:
            pdf_path.relative_to(expected_figure_root)
        except ValueError as error:
            raise ReleaseContractError(
                f"panel PDF is outside candidate misc root: {identity}"
            ) from error
        payload = build_pdf(row["panel_title"], source_rows)
        if (
            not payload.startswith(b"%PDF-1.4")
            or not payload.rstrip().endswith(b"%%EOF")
            or b"/BaseFont /Helvetica" not in payload
            or b"/F1 6 Tf" not in payload
        ):
            raise ReleaseContractError(f"rendered PDF contract failed: {identity}")
        prepared_panels.append((row, identity, pdf_path, payload))

    # All blueprints, source hashes, identities, and PDF payload contracts are
    # checked before creating the external candidate figure root.
    figure_root = ensure_candidate_output_root(project, "figures")
    if figure_root != expected_figure_root:
        raise ReleaseContractError("candidate figure-root identity drift")
    manifest_rows = []
    for row, identity, pdf_path, payload in prepared_panels:
        materialize_immutable(
            pdf_path, payload, figure_root, f"{identity} individual PDF"
        )
        pdf_bytes = pdf_path.read_bytes()
        if (
            not pdf_bytes.startswith(b"%PDF-1.4")
            or not pdf_bytes.rstrip().endswith(b"%%EOF")
            or b"/BaseFont /Helvetica" not in pdf_bytes
            or b"/F1 6 Tf" not in pdf_bytes
        ):
            raise ReleaseContractError(f"rendered PDF contract failed: {identity}")
        manifest_rows.append(
            {
                **row,
                "pdf_sha256": sha256_file(pdf_path),
                "pdf_bytes": pdf_path.stat().st_size,
                "status": "rendered_all_prespecified_rows",
            }
        )
    forbidden = [
        path
        for path in figure_root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".png", ".svg"}
    ]
    if forbidden:
        raise ReleaseContractError(
            f"candidate figures contain PNG/SVG files: {forbidden}"
        )
    figure_manifest = candidate / "manifests/figure_source_manifest.tsv"
    materialize_immutable_tsv(
        figure_manifest,
        manifest_rows,
        FIGURE_SOURCE_MANIFEST_FIELDS,
        candidate,
        "figure source manifest",
    )
    state = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "state": "REL03_PANELS_FROZEN",
        "input_snapshot_spec_sha256": base["spec_hash"],
        "rel02_state_sha256": sha256_file(candidate / "manifests/rel02_state.json"),
        "figure_source_manifest_sha256": sha256_file(figure_manifest),
        "n_panels": len(manifest_rows),
        "n_main_figures": 5,
        "n_supplementary_panels": len(expected_supplementary),
        "normal_font_pt": 6,
        "format": "PDF",
        "myojin_role": "supplement_only",
        "fixture_mode": fixture_mode,
        "canonical_promotion_status": "not_promoted",
    }
    state_path = candidate / "manifests/rel03_state.json"
    materialize_immutable_json(state_path, state, candidate, "REL-03 state")
    if rel02_state["figure_count"] != 5:
        raise ReleaseContractError("REL-02/REL-03 figure count disagreement")
    return {
        "status": "REL03_PANELS_FROZEN",
        "candidate_id": CANDIDATE_ID,
        "n_panels": len(manifest_rows),
        "n_main_figures": 5,
        "n_supplementary_panels": len(expected_supplementary),
        "myojin_role": "supplement_only",
        "figure_manifest_sha256": sha256_file(figure_manifest),
        "canonical_promotion_status": "not_promoted",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--fixture-mode", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        result = render_candidate_panels(args.project_root, args.fixture_mode)
    except ReleaseContractError as error:
        raise SystemExit(f"REL03_BLOCKED: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
