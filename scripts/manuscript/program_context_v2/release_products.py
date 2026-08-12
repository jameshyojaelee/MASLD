#!/usr/bin/env python3
"""Candidate-only REL-02--04 contracts, paths, and deterministic helpers."""

from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from release_common import (
    CANDIDATE_ID,
    DOWNSTREAM_FIELDS,
    SNAPSHOT_MANIFEST_FIELDS,
    ReleaseContractError,
    assert_no_symlinks,
    candidate_root,
    clean_relative_path,
    is_relative_to,
    parse_bool,
    parse_nonnegative_int,
    read_tsv_exact,
    require_sha256,
    sha256_file,
)
from fibrosis_candidate_contract import (
    TRUE_KLEINER_TRANSITIONS,
    validate_frozen_fibrosis_bundle,
)


# `load_figure_data.R` defines FIG_OUT/FIG_MISC. Candidate panels mirror that
# convention but remain below the isolated misc/candidates branch.
CANDIDATE_FIGURE_REL = Path("figures/misc/candidates") / CANDIDATE_ID
CANDIDATE_MANUSCRIPT_REL = Path("docs/manuscript/candidates") / CANDIDATE_ID

BASE_SELECTION_FIELDS = (
    "selection_id",
    "coordinator",
    "date",
    "artifact_id",
    "source_path",
    "source_sha256",
    "source_bytes",
    "snapshot_relpath",
    "artifact_role",
    "allowed_wording",
    "prohibited_wording",
)
BASE_SNAPSHOT_FIELDS = (
    "artifact_id",
    "artifact_role",
    "source_path_provenance",
    "source_sha256",
    "source_bytes",
    "snapshot_path",
    "snapshot_sha256",
    "snapshot_bytes",
    "allowed_wording",
    "prohibited_wording",
)
BASE_DOWNSTREAM_FIELDS = (
    "consumer_id",
    "artifact_id",
    "snapshot_path",
    "sha256",
    "bytes",
)

SOURCE_ROW_FIELDS = (
    "record_id",
    "group_id",
    "artifact_id",
    "claim_id",
    "section_id",
    "figure_id",
    "panel_id",
    "panel_title",
    "panel_role",
    "label",
    "category",
    "number_role",
    "plot_role",
    "value",
    "display_value",
    "numerator",
    "denominator",
    "unit",
    "biological_unit",
    "model_contrast",
    "effect_unit",
    "p_value",
    "q_value",
    "evidence_status",
    "negative_adequacy_criterion",
    "tested_universe",
    "source_dependence",
    "discovery_sources",
    "evaluation_sources",
    "reuse_detail",
    "independence_boundary",
    "allowed_wording",
    "prohibited_wording",
    "claim_text",
    "next_experiment",
    "manuscript_included",
    "plot_order",
    "is_control",
    "source_snapshot_path",
    "source_sha256",
)

# Continuous NMF axes are retained only as a supplementary decomposition.  The
# exact generated-row language below is the primary firewall; this predicate is
# a secondary defense for prose outside those generated rows.
NMF_DISCRETE_CLASS_CLAIM_PATTERNS = tuple(
    re.compile(pattern, re.IGNORECASE)
    for pattern in (
        r"\bsubtyp(?:e|es|ing)\b",
        r"\b(?:discrete|reproducible|stable|molecular|patient)\s+(?:patient\s+)?(?:cluster|clusters|class|classes|group|groups)\b",
        r"\bhard\s+(?:cluster|clusters|class|classes|partition|partitions)\b",
        r"\b(?:stratif(?:y|ies|ied)|segment(?:s|ed)?)\s+patients?\b",
        r"\bpatients?\s+(?:stratification|segmentation)\b",
        r"\b(?:separat(?:e|es|ed|ing)|divid(?:e|es|ed|ing))\s+patients?\s+into\s+(?:stable\s+|molecular\s+)?(?:groups|clusters|classes)\b",
        r"\b(?:assigned|identified|recovered|defined|found|produced|yielded)\s+(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten)\s+(?:patient\s+)?(?:groups?|clusters?|classes?|subtypes?|partitions?)\b",
        r"\b(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten)\s+(?:stable\s+|molecular\s+|patient\s+)*(?:groups?|clusters?|classes?|subtypes?|partitions?)\b",
    )
)

NMF_S2A_ALLOWED_WORDING = (
    "continuous interpretive NMF axis with descriptive loading distribution"
)
NMF_S2A_PROHIBITED_WORDING = (
    "subtype, hard cluster, patient class, reproducible stratification, "
    "transition, or biomarker"
)
NMF_S2B_ALLOWED_WORDING = "factor-axis stability with continuous interpretation only"
NMF_S2B_PROHIBITED_WORDING = (
    "stable subtype, reproducible hard partition, patient class, clinical "
    "stratifier, or biomarker"
)
NMF_NEXT_EXPERIMENT = (
    "prospective outcome-linked validation before any patient-stratification use"
)
NMF_S2A_RECORD_RE = re.compile(r"^nmf_k(?P<k>[46])_p(?P<program>[1-6])_continuous$")
NMF_S2B_RECORD_RE = re.compile(
    r"^nmf_k(?P<k>[46])_(?P<metric>mean_cophenetic|mean_silhouette|"
    r"mean_matched_cosine|min_top50_jaccard)$"
)


def contains_nmf_discrete_class_claim(text: str) -> bool:
    """Return True when NMF prose implies reproducible hard patient classes."""

    for pattern in NMF_DISCRETE_CLASS_CLAIM_PATTERNS:
        for match in pattern.finditer(text):
            # Boundary statements such as "not used as a patient
            # stratification system" are required and must not be confused
            # with a positive subtype claim.  Negation must precede the matched
            # phrase in the same short clause; distant/global disclaimers do
            # not immunize a positive claim elsewhere in the row.
            prefix = text[max(0, match.start() - 64) : match.start()]
            suffix = text[match.end() : match.end() + 64]
            preceding_boundary = re.search(
                r"\b(?:no|not|never|cannot|can't|do not|does not|is not|are not|"
                r"prohibit(?:ed)?|avoid(?:ed)?)\b[^,.;:]{0,48}$",
                prefix,
                re.IGNORECASE,
            )
            following_boundary = re.match(
                r"^[^,.;:]{0,32}\b(?:is|are|was|were|do|does)\s+not\s+"
                r"(?:inferred|assigned|defined|established|supported|claimed)\b",
                suffix,
                re.IGNORECASE,
            )
            numeric_assignment = re.search(
                r"\b(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten)\b",
                match.group(0),
                re.IGNORECASE,
            )
            if (preceding_boundary or following_boundary) and not numeric_assignment:
                continue
            return True
    return False


def validate_nmf_source_language(row: Mapping[str, str]) -> None:
    """Require one of the two exact approved generated-row templates."""

    record_id = row["record_id"]
    if row["panel_id"] == "S2A":
        match = NMF_S2A_RECORD_RE.fullmatch(record_id)
        if match is None or int(match.group("program")) > int(match.group("k")):
            raise ReleaseContractError(
                f"{record_id} is not an approved continuous-axis record ID"
            )
        expected = {
            "number_role": "median_continuous_loading",
            "model_contrast": "descriptive sample-level continuous loading distribution",
            "effect_unit": "median continuous loading; not an inferential effect or class assignment",
            "allowed_wording": NMF_S2A_ALLOWED_WORDING,
            "prohibited_wording": NMF_S2A_PROHIBITED_WORDING,
            "claim_text": (
                f"k={match.group('k')} P{match.group('program')} is retained only as a "
                "continuous interpretive axis across the frozen sample set."
            ),
            "next_experiment": NMF_NEXT_EXPERIMENT,
        }
    elif row["panel_id"] == "S2B":
        match = NMF_S2B_RECORD_RE.fullmatch(record_id)
        if match is None:
            raise ReleaseContractError(
                f"{record_id} is not an approved factor-stability record ID"
            )
        expected = {
            "number_role": match.group("metric"),
            "model_contrast": "three-seed internal stability audit",
            "effect_unit": "factor-stability statistic; not a sample-partition reproducibility estimate",
            "allowed_wording": NMF_S2B_ALLOWED_WORDING,
            "prohibited_wording": NMF_S2B_PROHIBITED_WORDING,
            "claim_text": (
                f"The k={match.group('k')} factor audit is shown with an explicit "
                "boundary that factor stability does not imply stable patient partitions."
            ),
            "next_experiment": NMF_NEXT_EXPERIMENT,
        }
    else:
        raise ReleaseContractError(
            f"{record_id} uses an unapproved NMF supplementary panel"
        )
    drift = {
        field: {"expected": value, "observed": row.get(field)}
        for field, value in expected.items()
        if row.get(field) != value
    }
    if drift:
        raise ReleaseContractError(
            f"{record_id} differs from the approved NMF language template: {drift}"
        )


NUMBERS_FIELDS = (
    "number_id",
    "claim_id",
    "number_role",
    "raw_value",
    "display_value",
    "numerator",
    "denominator",
    "unit_of_replication",
    "model_contrast",
    "source_artifact",
    "source_row_filter",
    "release_id",
    "allowed_wording",
    "prohibited_wording",
    "figure_id",
    "panel_id",
    "manuscript_included",
)
CLAIM_FIELDS = (
    "claim_id",
    "section_id",
    "claim_text",
    "claim_status",
    "number_ids",
    "source_artifact_ids",
    "allowed_wording",
    "prohibited_wording",
    "source_dependence",
    "tested_universe",
    "biological_unit",
)
SOURCE_DEPENDENCY_FIELDS = (
    "claim_id",
    "source_dependence_class",
    "discovery_source",
    "evaluation_source",
    "reuse_detail",
    "independence_boundary",
)
SOURCE_DEPENDENCE_CLASSES = {
    "independent",
    "reused_source",
    "partially_dependent",
    "source_dependent",
    "mixed",
}
SOURCE_ID_RE = re.compile(r"^[A-Za-z0-9_.:-]+$")
ACCEPTANCE_GATE_FIELDS = (
    "gate_id",
    "workstream_id",
    "terminal_state",
    "gate_verdict",
    "include_main",
    "include_supplement",
    "figure_role",
    "closure_sha256",
)
ANALYSIS_MANIFEST_FIELDS = (
    "artifact_id",
    "artifact_role",
    "candidate_path",
    "sha256",
    "bytes",
    "producer_id",
    "input_snapshot_spec_sha256",
)
PANEL_BLUEPRINT_FIELDS = (
    "figure_id",
    "panel_id",
    "figure_order",
    "panel_order",
    "panel_title",
    "panel_role",
    "source_table_path",
    "source_table_sha256",
    "pdf_path",
    "normal_font_pt",
    "format",
)
FIGURE_SOURCE_MANIFEST_FIELDS = PANEL_BLUEPRINT_FIELDS + (
    "pdf_sha256",
    "pdf_bytes",
    "status",
)
MANUSCRIPT_MANIFEST_FIELDS = (
    "document_id",
    "candidate_path",
    "sha256",
    "bytes",
    "journal_branch",
    "figure_count",
    "myojin_role",
    "number_ids",
    "claim_ids",
)
TRANSITION_PRODUCT_FIELDS = (
    "product_id",
    "phase",
    "artifact_role",
    "project_relative_path",
    "sha256",
    "bytes",
    "producer_id",
    "input_snapshot_spec_sha256",
)

EXPECTED_FIGURES = ("Figure1", "Figure2", "Figure3", "Figure4", "Figure5")
EXPECTED_FIGURE_TITLES = {
    "Figure1": "Dataset interface and evidence observability",
    "Figure2": "Sample- and donor-resolved established-state transcriptomics",
    "Figure3": "Genetics and regulatory context",
    "Figure4": "Physical context, assay observability, and prespecified external challenges",
    "Figure5": "MASLD Gene Catalog and translational boundaries",
}
EXPECTED_SUPPLEMENTARY_FIGURES = ("FigureS1", "FigureS2")
EXPECTED_SUPPLEMENTARY_FIGURE_TITLES = {
    "FigureS1": "Prespecified Myojin HLF functional stress test",
    "FigureS2": "Continuous NMF axes and factor stability",
}
EXPECTED_SUPPLEMENTARY_PANELS = {
    "FigureS1": (
        (
            "S1A",
            "Complete Myojin evidence-class stress test",
            "supplementary_functional",
        ),
        (
            "S1B",
            "Complete Myojin program stress test",
            "supplementary_functional",
        ),
    ),
    "FigureS2": (
        (
            "S2A",
            "Continuous k4/k6 NMF loading distributions",
            "supplementary_continuous_programs",
        ),
        (
            "S2B",
            "Three-seed factor stability and hard-partition boundary",
            "supplementary_continuous_programs",
        ),
    ),
}
EXPECTED_RESULTS_SECTIONS = (
    "resource_interface",
    "established_state_transcriptomics",
    "genetics_context",
    "physical_context",
    "evidence_passports",
)
EXPECTED_SUPPLEMENTARY_SECTIONS = (
    "supplementary_myojin",
    "supplementary_nmf",
)
JOURNAL_BRANCH = "Cell Genomics Resource"
MYOJIN_ROLE = "supplement_only"
TRUE_KLEINER_TRANSITION_CENSUS = {
    contrast: (int(spec["n_samples"]), len(spec["cohorts"]))
    for contrast, spec in TRUE_KLEINER_TRANSITIONS.items()
}


def candidate_figure_root(project_root: Path) -> Path:
    return project_root.resolve() / CANDIDATE_FIGURE_REL


def candidate_manuscript_root(project_root: Path) -> Path:
    return project_root.resolve() / CANDIDATE_MANUSCRIPT_REL


def assert_exact_output_root(project_root: Path, proposed: Path, scope: str) -> Path:
    project = project_root.resolve()
    expected_by_scope = {
        "release": candidate_root(project),
        "figures": candidate_figure_root(project),
        "manuscript": candidate_manuscript_root(project),
    }
    if scope not in expected_by_scope:
        raise ReleaseContractError(f"unknown candidate output scope: {scope}")
    expected = Path(os.path.abspath(expected_by_scope[scope]))
    observed = Path(os.path.abspath(proposed))
    if observed != expected:
        raise ReleaseContractError(
            f"{scope} output root must be exactly {expected}; observed {observed}"
        )
    cursor = project
    for part in observed.relative_to(project).parts:
        cursor /= part
        if cursor.is_symlink():
            raise ReleaseContractError(
                f"{scope} output root uses a symlink component: {cursor}"
            )
    return observed


def ensure_candidate_output_root(project_root: Path, scope: str) -> Path:
    root = assert_exact_output_root(
        project_root,
        {
            "release": candidate_root(project_root),
            "figures": candidate_figure_root(project_root),
            "manuscript": candidate_manuscript_root(project_root),
        }[scope],
        scope,
    )
    if root.exists() and (not root.is_dir() or root.is_symlink()):
        raise ReleaseContractError(f"candidate {scope} root is unsafe: {root}")
    root.mkdir(parents=True, exist_ok=True)
    assert_no_symlinks(root)
    return root


def read_json_object(path: Path, context: str) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(f"{context} is missing or symlinked: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ReleaseContractError(f"{context} is invalid JSON: {path}") from error
    if not isinstance(payload, dict):
        raise ReleaseContractError(f"{context} must be a JSON object: {path}")
    return payload


def parse_float(value: str, context: str, allow_blank: bool = False) -> float | None:
    if allow_blank and not value.strip():
        return None
    try:
        result = float(value)
    except ValueError as error:
        raise ReleaseContractError(f"{context} is not numeric: {value!r}") from error
    if result != result or result in {float("inf"), float("-inf")}:
        raise ReleaseContractError(f"{context} is non-finite: {value!r}")
    return result


def parse_source_ids(value: str, context: str) -> set[str]:
    """Parse a deterministic, nonempty semicolon-delimited source-ID set."""

    parts = value.split(";") if value else []
    if (
        not parts
        or any(not SOURCE_ID_RE.fullmatch(part) for part in parts)
        or parts != sorted(set(parts))
    ):
        raise ReleaseContractError(
            f"{context} must be a sorted, unique, semicolon-delimited set of "
            f"stable source IDs: {value!r}"
        )
    return set(parts)


def validate_source_dependency_contract(
    dependence_class: str,
    discovery_text: str,
    evaluation_text: str,
    reuse_detail: str,
    independence_boundary: str,
    context: str,
) -> None:
    """Reject independence claims that are inconsistent with named source sets."""

    if dependence_class not in SOURCE_DEPENDENCE_CLASSES:
        raise ReleaseContractError(
            f"{context} has unsupported source-dependence class: {dependence_class!r}"
        )
    discovery = parse_source_ids(discovery_text, f"{context} discovery_sources")
    evaluation = parse_source_ids(evaluation_text, f"{context} evaluation_sources")
    if not reuse_detail.strip() or not independence_boundary.strip():
        raise ReleaseContractError(
            f"{context} requires nonblank reuse_detail and independence_boundary"
        )
    overlap = discovery & evaluation
    nonoverlap = discovery ^ evaluation
    if dependence_class == "independent" and overlap:
        raise ReleaseContractError(
            f"{context} labels overlapping sources independent: {sorted(overlap)}"
        )
    if dependence_class in {"reused_source", "source_dependent"} and not overlap:
        raise ReleaseContractError(
            f"{context} declares source reuse without a named overlapping source"
        )
    if dependence_class == "partially_dependent" and (not overlap or not nonoverlap):
        raise ReleaseContractError(
            f"{context} partially_dependent requires both shared and distinct sources"
        )


def read_source_rows(path: Path) -> list[dict[str, str]]:
    rows = read_tsv_exact(path, SOURCE_ROW_FIELDS[:-2])
    if not rows:
        raise ReleaseContractError(f"candidate evidence source is empty: {path}")
    seen = set()
    for row in rows:
        record_id = row["record_id"].strip()
        if not record_id or record_id in seen:
            raise ReleaseContractError(
                f"blank or duplicate source record_id in {path}: {record_id!r}"
            )
        seen.add(record_id)
        for field in (
            "group_id",
            "artifact_id",
            "claim_id",
            "section_id",
            "panel_id",
            "label",
            "number_role",
            "plot_role",
            "display_value",
            "unit",
            "biological_unit",
            "model_contrast",
            "evidence_status",
            "tested_universe",
            "source_dependence",
            "discovery_sources",
            "evaluation_sources",
            "reuse_detail",
            "independence_boundary",
            "allowed_wording",
            "prohibited_wording",
            "claim_text",
        ):
            if not row[field].strip():
                raise ReleaseContractError(
                    f"source row {record_id} has blank required field {field}"
                )
        validate_source_dependency_contract(
            row["source_dependence"],
            row["discovery_sources"],
            row["evaluation_sources"],
            row["reuse_detail"],
            row["independence_boundary"],
            f"source row {record_id}",
        )
        parse_float(row["value"], f"{record_id} value")
        parse_float(row["p_value"], f"{record_id} p_value", allow_blank=True)
        parse_float(row["q_value"], f"{record_id} q_value", allow_blank=True)
        parse_nonnegative_int(row["plot_order"], f"{record_id} plot_order")
        parse_bool(row["manuscript_included"], f"{record_id} manuscript_included")
        parse_bool(row["is_control"], f"{record_id} is_control")
        if row["evidence_status"] == "tested_negative":
            if not row["negative_adequacy_criterion"].strip() or not any(
                token in row["allowed_wording"].lower()
                for token in ("adequate-negative", "equivalence criterion")
            ):
                raise ReleaseContractError(
                    f"{record_id} uses tested_negative without an explicit "
                    "adequate-negative/equivalence criterion"
                )
        elif row["negative_adequacy_criterion"].strip():
            raise ReleaseContractError(
                f"{record_id} carries a negative-adequacy criterion but is not tested_negative"
            )
        combined_text = " ".join(
            row[field]
            for field in (
                "label",
                "category",
                "number_role",
                "model_contrast",
                "effect_unit",
                "allowed_wording",
                "claim_text",
            )
        ).lower()
        if "broad validation" in combined_text:
            raise ReleaseContractError(
                f"{record_id} overstates an external challenge as broad validation"
            )
        if any(token in combined_text for token in ("vu array", "cosmx array")):
            if (
                "donor unknown" not in row["biological_unit"].lower()
                or "technical" not in row["unit"].lower()
            ):
                raise ReleaseContractError(
                    f"{record_id} must keep unresolved arrays as technical units with donor n unknown"
                )
        effect_text = row["effect_unit"].lower()
        if "moran" in combined_text and any(
            token in effect_text and f"not {token}" not in effect_text
            for token in ("sampling se", "sampling ci")
        ):
            raise ReleaseContractError(
                f"{record_id} mislabels matched-null Moran dispersion as sampling uncertainty"
            )
        if row["section_id"] == "supplementary_myojin":
            if row["figure_id"] != "FigureS1":
                raise ReleaseContractError(
                    f"{record_id} places Myojin evidence outside FigureS1"
                )
            if any(token in combined_text for token in ("validated", "precise null")):
                raise ReleaseContractError(
                    f"{record_id} overstates the Myojin assay-specific result"
                )
            if not any(
                token in combined_text
                for token in ("assay-specific non-support", "nonconfirmatory")
            ):
                raise ReleaseContractError(
                    f"{record_id} lacks assay-specific non-support/nonconfirmatory wording"
                )
        if row["section_id"] == "supplementary_nmf":
            if row["figure_id"] != "FigureS2":
                raise ReleaseContractError(
                    f"{record_id} places continuous NMF evidence outside FigureS2"
                )
            validate_nmf_source_language(row)
            if contains_nmf_discrete_class_claim(combined_text):
                raise ReleaseContractError(
                    f"{record_id} overstates continuous NMF axes as subtypes"
                )
        if row["figure_id"]:
            if row["figure_id"] not in EXPECTED_FIGURES and not row[
                "figure_id"
            ].startswith("FigureS"):
                raise ReleaseContractError(
                    f"source row {record_id} has invalid figure_id {row['figure_id']}"
                )
        if row["section_id"] not in (
            *EXPECTED_RESULTS_SECTIONS,
            *EXPECTED_SUPPLEMENTARY_SECTIONS,
        ):
            raise ReleaseContractError(
                f"source row {record_id} has invalid section_id {row['section_id']}"
            )
    return rows


def write_tsv_payload(
    rows: Iterable[Mapping[str, object]], fields: Sequence[str]
) -> bytes:
    from rel01_snapshot_candidate import canonical_tsv_bytes

    return canonical_tsv_bytes(list(rows), tuple(fields))


def project_relative_str(project_root: Path, path: Path) -> str:
    project = project_root.resolve()
    resolved = path.resolve()
    if not is_relative_to(resolved, project):
        raise ReleaseContractError(f"product escapes project root: {path}")
    return resolved.relative_to(project).as_posix()


def verify_file_row(
    root: Path,
    relative_text: str,
    expected_sha256: str,
    expected_bytes: str,
    context: str,
) -> Path:
    relative = clean_relative_path(relative_text, context)
    path = root / relative
    require_sha256(expected_sha256, f"{context} sha256")
    size = parse_nonnegative_int(expected_bytes, f"{context} bytes")
    if (
        not path.is_file()
        or path.is_symlink()
        or path.stat().st_size != size
        or sha256_file(path) != expected_sha256
    ):
        raise ReleaseContractError(f"{context} hash/byte mismatch: {path}")
    return path


def read_table_flexible(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(f"input table is missing or symlinked: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(
            handle, delimiter="\t" if path.suffix == ".tsv" else ","
        )
        fields = list(reader.fieldnames or [])
        return fields, [dict(row) for row in reader]


def validate_base_extension(
    project_root: Path,
    candidate: Path,
    spec_hash: str,
    fixture_mode: bool,
) -> dict[str, dict[str, str]]:
    transition_path = candidate / "manifests/base_input_transition.json"
    transition = read_json_object(transition_path, "base input transition")
    expected = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "transition": "BASE_INPUTS_FROZEN",
        "input_snapshot_spec_sha256": spec_hash,
        "fixture_mode": fixture_mode,
        "scientific_assembly_performed": False,
        "canonical_promotion_status": "not_promoted",
    }
    for key, value in expected.items():
        if transition.get(key) != value:
            raise ReleaseContractError(
                f"base input transition drift for {key}: {transition.get(key)!r}"
            )
    state_path = candidate / "manifests/candidate_state.json"
    if transition.get("from_candidate_state_sha256") != sha256_file(state_path):
        raise ReleaseContractError(
            "base input transition does not preserve immutable candidate state"
        )
    selection_path = candidate / "manifests/base_input_selection.tsv"
    selection_rows = read_tsv_exact(selection_path, BASE_SELECTION_FIELDS)
    from rebuild_source_spec import load_rebuild_source_spec

    rebuild_spec, _ = load_rebuild_source_spec(candidate, fixture_mode=fixture_mode)
    normalized_selection = sorted(selection_rows, key=lambda row: row["artifact_id"])
    if transition.get("selection_sha256") != rebuild_spec.get(
        "base_selection_sha256"
    ) or normalized_selection != rebuild_spec.get("base_selection"):
        raise ReleaseContractError(
            "BASE transition differs from the frozen rebuild source identity"
        )
    if transition.get("selection_sha256") != sha256_file(selection_path):
        raise ReleaseContractError("base input selection hash drift")
    if transition.get("selection_snapshot_sha256") != sha256_file(selection_path):
        raise ReleaseContractError("copied base selection hash drift")
    snapshot_path = candidate / "manifests/base_input_snapshot_manifest.tsv"
    downstream_path = candidate / "manifests/base_downstream_inputs.tsv"
    if transition.get("base_snapshot_manifest_sha256") != sha256_file(snapshot_path):
        raise ReleaseContractError("base snapshot manifest hash drift")
    if transition.get("base_downstream_manifest_sha256") != sha256_file(
        downstream_path
    ):
        raise ReleaseContractError("base downstream manifest hash drift")
    snapshot_rows = read_tsv_exact(snapshot_path, BASE_SNAPSHOT_FIELDS)
    downstream_rows = read_tsv_exact(downstream_path, BASE_DOWNSTREAM_FIELDS)
    if transition.get("n_base_inputs") != len(snapshot_rows):
        raise ReleaseContractError("base transition input count drift")
    selection_by_id = {row["artifact_id"]: row for row in selection_rows}
    snapshot_by_id = {row["artifact_id"]: row for row in snapshot_rows}
    downstream_by_id = {row["artifact_id"]: row for row in downstream_rows}
    if (
        len(selection_by_id) != len(selection_rows)
        or len(snapshot_by_id) != len(snapshot_rows)
        or len(downstream_by_id) != len(downstream_rows)
        or set(selection_by_id) != set(snapshot_by_id)
        or set(selection_by_id) != set(downstream_by_id)
    ):
        raise ReleaseContractError(
            "base selection, snapshot, and downstream artifact identities differ"
        )
    observed_paths = set()
    index: dict[str, dict[str, str]] = {}
    for artifact_id, selected in selection_by_id.items():
        snap = snapshot_by_id[artifact_id]
        downstream = downstream_by_id[artifact_id]
        source_provenance = clean_relative_path(
            snap["source_path_provenance"],
            f"base {artifact_id} source provenance",
        )
        candidate_relative = project_relative_str(project_root, candidate)
        if source_provenance == candidate_relative or source_provenance.startswith(
            f"{candidate_relative}/"
        ):
            raise ReleaseContractError(
                f"base source provenance points into candidate: {artifact_id}"
            )
        relative = clean_relative_path(
            snap["snapshot_path"], f"base {artifact_id} snapshot path"
        )
        if not relative.startswith("inputs/BASE/"):
            raise ReleaseContractError(
                f"base snapshot escaped inputs/BASE/: {artifact_id}"
            )
        destination = candidate / relative
        source_hash = require_sha256(
            snap["source_sha256"], f"base {artifact_id} source hash"
        )
        source_bytes = parse_nonnegative_int(
            snap["source_bytes"], f"base {artifact_id} source bytes"
        )
        if (
            not destination.is_file()
            or destination.is_symlink()
            or destination.stat().st_size != source_bytes
            or sha256_file(destination) != source_hash
        ):
            raise ReleaseContractError(f"base copied input drift: {artifact_id}")
        expected_pairs = {
            "artifact_role": selected["artifact_role"],
            "source_sha256": selected["source_sha256"],
            "source_bytes": selected["source_bytes"],
            "snapshot_path": selected["snapshot_relpath"],
            "snapshot_sha256": selected["source_sha256"],
            "snapshot_bytes": selected["source_bytes"],
            "allowed_wording": selected["allowed_wording"],
            "prohibited_wording": selected["prohibited_wording"],
            "source_path_provenance": selected["source_path"],
        }
        for field, expected_value in expected_pairs.items():
            if snap[field] != expected_value:
                raise ReleaseContractError(
                    f"base snapshot contract drift for {artifact_id} {field}"
                )
        if (
            downstream["consumer_id"] != f"BASE:{artifact_id}"
            or downstream["snapshot_path"] != relative
            or downstream["sha256"] != source_hash
            or downstream["bytes"] != str(source_bytes)
        ):
            raise ReleaseContractError(f"base downstream contract drift: {artifact_id}")
        observed_paths.add(relative)
        index[f"BASE:{artifact_id}"] = {
            "source_key": f"BASE:{artifact_id}",
            "workstream_id": "BASE",
            "artifact_id": artifact_id,
            "artifact_role": snap["artifact_role"],
            "snapshot_path": relative,
            "sha256": source_hash,
            "bytes": str(source_bytes),
            "allowed_wording": snap["allowed_wording"],
            "prohibited_wording": snap["prohibited_wording"],
        }
    tree_paths = {
        path.relative_to(candidate).as_posix()
        for path in (candidate / "inputs/BASE").rglob("*")
        if path.is_file()
    }
    if tree_paths != observed_paths:
        raise ReleaseContractError(
            "inputs/BASE tree does not exactly match its frozen manifest"
        )
    validate_frozen_fibrosis_bundle(candidate, fixture_mode=fixture_mode)
    return index


def load_frozen_input_index(
    project_root: Path,
    candidate: Path,
    spec_hash: str,
    fixture_mode: bool,
) -> dict[str, dict[str, str]]:
    snapshot_rows = read_tsv_exact(
        candidate / "manifests/input_snapshot_manifest.tsv",
        SNAPSHOT_MANIFEST_FIELDS,
    )
    snapshot_by_path = {row["snapshot_relpath"]: row for row in snapshot_rows}
    downstream_rows = read_tsv_exact(
        candidate / "manifests/downstream_inputs.tsv", DOWNSTREAM_FIELDS
    )
    index: dict[str, dict[str, str]] = {}
    for row in downstream_rows:
        relative = clean_relative_path(
            row["snapshot_path"], f"{row['consumer_id']} snapshot path"
        )
        snapshot = snapshot_by_path.get(relative)
        if snapshot is None:
            raise ReleaseContractError(
                f"downstream input lacks snapshot record: {relative}"
            )
        key = f"{row['workstream_id']}:{row['artifact_id']}"
        if key in index:
            raise ReleaseContractError(f"duplicate frozen source key: {key}")
        index[key] = {
            "source_key": key,
            "workstream_id": row["workstream_id"],
            "artifact_id": row["artifact_id"],
            "artifact_role": snapshot["artifact_role"],
            "snapshot_path": relative,
            "sha256": row["sha256"],
            "bytes": row["bytes"],
            "allowed_wording": "",
            "prohibited_wording": "",
        }
    base_index = validate_base_extension(
        project_root, candidate, spec_hash, fixture_mode
    )
    overlap = set(index) & set(base_index)
    if overlap:
        raise ReleaseContractError(f"base/workstream source-key collision: {overlap}")
    index.update(base_index)
    return index
