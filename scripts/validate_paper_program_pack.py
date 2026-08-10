#!/usr/bin/env python3
"""Validate the authoritative 2026-08-07 paper-program documentation pack."""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path


PACK_NAMES = (
    "00_INDEX_AND_ORCHESTRATION.md",
    "10_SPATIAL_ACQUISITION_AND_SOURCE_GATES.md",
    "11_YAKUBOVSKY_EARLY_LIPID_CONTEXT.md",
    "12_GSE287826_DONOR_REPLICATION.md",
    "13_SPATIAL_CONTEXT_AND_FIG4_INTEGRATION.md",
    "20_HOTSPOT_REGISTRY_AND_FIG2.md",
    "30_GENETICS_CONTEXT_AND_FIG3.md",
    "40_MYOJIN_FUNCTIONAL_STRESS_TEST.md",
    "41_PUBLIC_FUNCTIONAL_INDUCTION_AND_REVERSAL.md",
    "42_RISK_TO_STATE_FUNCTIONAL_DOUBLE_DISSOCIATION.md",
    "43_CHRONIC_STATE_REVERSAL_AND_REGULATORY_RISK_BRIDGE.md",
    "44_MULTICELLULAR_STATE_ASSEMBLY_AND_RESPONSE_TOPOLOGY.md",
    "45_SOURCE_INDEPENDENT_RISK_TO_STATE_RELAY.md",
    "45A_EXPERIMENTAL_PLATFORM_ROUTING_AND_OUTREACH.md",
    "46_FULL_GOAL_COMPLETION_NEW_DATA_PROTOCOL.md",
    "46A_CLINICAL_COLLABORATOR_ROUTING_AND_OUTREACH.md",
    "46B_CLINICAL_RESPONSE_INTAKE_AND_TWO_COHORT_GATE.md",
    "46C_BLINDED_CLINICAL_PREFLIGHT_AND_MODEL_FREEZE.md",
    "50_EVIDENCE_PASSPORTS_AND_PORTAL.md",
    "60_FIGURE_MANUSCRIPT_AND_RELEASE.md",
)
HEADER_KEYS = (
    "status",
    "authority",
    "owner",
    "upstream_inputs",
    "downstream_consumers",
    "start_gate",
    "completion_gate",
    "owned_output_paths",
    "forbidden_output_paths",
    "supersedes",
    "last_decision_date",
)
PREDECESSORS = (
    "2026-08-07_competitor_response_computational_article_pivot.md",
    "2026-08-07_competitor_response_final_direction.md",
    "2026-08-07_competitor_response_reconciled_plan.md",
)
HEADER_RE = re.compile(r"^([a-z][a-z0-9_]*):(?:\s*(.*))?$")
BACKTICK_MD_RE = re.compile(r"`([^`\n]+\.md)`")
MARKDOWN_LINK_RE = re.compile(r"\[[^\]]*\]\(([^)]+\.md)\)")


class ValidationError(RuntimeError):
    """One or more plan-pack invariants failed."""


def extract_header(text: str, path: Path) -> tuple[list[str], dict[str, list[str]]]:
    match = re.search(r"```yaml\n(.*?)\n```", text, flags=re.DOTALL)
    if match is None:
        raise ValidationError(f"{path}: missing fenced YAML execution header")
    keys: list[str] = []
    values: dict[str, list[str]] = defaultdict(list)
    current: str | None = None
    for line in match.group(1).splitlines():
        top = HEADER_RE.match(line)
        if top:
            current = top.group(1)
            keys.append(current)
            inline = (top.group(2) or "").strip()
            if inline:
                values[current].append(inline)
        elif current is not None and line.startswith("  - "):
            values[current].append(line[4:].strip())
        elif line.strip() and current is not None:
            values[current].append(line.strip())
    return keys, values


def resolve_reference(root: Path, pack: Path, source: Path, reference: str) -> bool:
    reference = reference.strip().strip("<>")
    if "#" in reference:
        reference = reference.split("#", 1)[0]
    candidates = [
        source.parent / reference,
        root / reference,
        root / "docs" / reference,
        root / "docs" / "manuscript" / reference,
        pack / Path(reference).name,
    ]
    return any(candidate.resolve().is_file() for candidate in candidates)


def assert_acyclic(edges: dict[str, set[str]]) -> None:
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(node: str, trail: list[str]) -> None:
        if node in visiting:
            cycle = " -> ".join([*trail, node])
            raise ValidationError(f"workstream dependency cycle: {cycle}")
        if node in visited:
            return
        visiting.add(node)
        for downstream in sorted(edges.get(node, set())):
            visit(downstream, [*trail, node])
        visiting.remove(node)
        visited.add(node)

    for node in PACK_NAMES:
        visit(node, [])


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    pack = root / "docs" / "plans" / "2026-08-07_paper_program"
    observed = tuple(sorted(path.name for path in pack.glob("*.md")))
    expected = tuple(sorted(PACK_NAMES))
    errors: list[str] = []
    if observed != expected:
        errors.append(
            f"pack file universe drift: observed={observed!r}; expected={expected!r}"
        )

    headers: dict[str, dict[str, list[str]]] = {}
    owner_paths: dict[str, list[str]] = defaultdict(list)
    edges: dict[str, set[str]] = defaultdict(set)
    for name in PACK_NAMES:
        path = pack / name
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        try:
            keys, values = extract_header(text, path)
        except ValidationError as error:
            errors.append(str(error))
            continue
        headers[name] = values
        if tuple(keys) != HEADER_KEYS:
            errors.append(f"{name}: header key/order drift: {tuple(keys)!r}")
        decision_values = values.get("last_decision_date", [])
        try:
            decision_date = date.fromisoformat(decision_values[0]) if len(decision_values) == 1 else None
        except ValueError:
            decision_date = None
        if decision_date is None or decision_date < date(2026, 8, 7) or decision_date > date.today():
            errors.append(f"{name}: invalid last_decision_date {decision_values!r}")
        for owned in values.get("owned_output_paths", []):
            owner_paths[owned].append(name)
        for upstream in values.get("upstream_inputs", []):
            upstream_name = Path(upstream.strip("` ")).name
            if upstream_name in PACK_NAMES:
                edges[upstream_name].add(name)

        references = set(MARKDOWN_LINK_RE.findall(text))
        references.update(
            reference
            for reference in BACKTICK_MD_RE.findall(text)
            if "/" in reference or Path(reference).name in PACK_NAMES
        )
        for reference in sorted(references):
            if not resolve_reference(root, pack, path, reference):
                errors.append(f"{name}: unresolved Markdown reference {reference!r}")

    duplicate_owners = {
        output: owners for output, owners in owner_paths.items() if len(owners) != 1
    }
    if duplicate_owners:
        errors.append(f"exact owned-output collisions: {duplicate_owners!r}")
    try:
        assert_acyclic(edges)
    except ValidationError as error:
        errors.append(str(error))

    for name in PREDECESSORS:
        path = root / "docs" / "plans" / name
        if not path.is_file():
            errors.append(f"missing retained predecessor: {path}")
            continue
        prefix = "\n".join(path.read_text(encoding="utf-8").splitlines()[:10]).lower()
        if "authority notice" not in prefix or "superseded for execution" not in prefix:
            errors.append(f"{name}: missing top-of-file supersession authority banner")

    if errors:
        print("PLAN_PACK_VALIDATION_FAIL", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        raise SystemExit(1)
    print(
        "PLAN_PACK_VALIDATION_PASS\t"
        f"documents={len(PACK_NAMES)}\t"
        f"owned_paths={len(owner_paths)}\t"
        f"dependency_edges={sum(map(len, edges.values()))}\t"
        f"predecessor_banners={len(PREDECESSORS)}"
    )


if __name__ == "__main__":
    main()
