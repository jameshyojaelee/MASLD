#!/usr/bin/env python3
"""Record the zero-applicable-conditions rule as a CAMPAIGN-WIDE check convention.

W3's frozen record states the finding, but states it as something preserved for
that lane. It is not a W3 rule. It is a property of every check that applies
not_applicable honestly, and it is the fifth instance of one failure shape.

The sharpest part is the part the W3 record does not carry: this trap was
CREATED by the fix for the earlier ones. Before honest not_applicable, a check
in this state reported a misleading met/total. After, it can report a spurious
PASS. The correctness fix moved the failure rather than removing it, and moved
it somewhere quieter - every individual condition behaves correctly and the bug
lives only in the aggregation.

Additive. Modifies nothing.
"""

from __future__ import annotations

import argparse, hashlib, json
from pathlib import Path


def sha256_file(p: Path) -> str:
    d = hashlib.sha256()
    with p.open("rb") as h:
        for b in iter(lambda: h.read(8 << 20), b""):
            d.update(b)
    return d.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--w3-record", required=True, type=Path)
    ap.add_argument("--w2-gate", required=True, type=Path)
    ap.add_argument("--typed-evidence-graph", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    from masld_bench.artifacts import freeze_tree

    if a.output.exists():
        raise RuntimeError("refusing to overwrite")

    rec = {
        "schema_version": "masld-bench-gate-convention-v1",
        "record_id": "zero_applicable_conditions_gate_convention_v1",
        "status": "campaign_wide_convention_applies_to_every_future_gate",
        "scope": "CAMPAIGN-WIDE. This is not a W3 rule and must not be read as one.",
        "record_is_additive_overlay": True,
        "modifies_nothing": True,
        "amends": {
            "w3_record_artifacts_sha256": sha256_file(a.w3_record / "ARTIFACTS.json"),
            "reason": (
                "The W3 record states the finding but frames it as preserved for that lane. "
                "The rule is a property of every gate, so it is restated here at campaign scope."
            ),
        },

        "the_rule": {
            "statement": (
                "Every evaluator must handle zero applicable conditions EXPLICITLY, and must "
                "never let a boolean aggregation over an empty set stand for a verdict."
            ),
            "the_mechanism": "all([]) returns True in Python, and the equivalent holds in every language with a fold over an empty conjunction.",
            "the_consequence": (
                "An evaluator that applies not_applicable honestly, reaches zero applicable "
                "conditions, and then computes all(c.met for c in applicable) will emit PASS "
                "having tested nothing."
            ),
            "required_third_state": "NO_APPLICABLE_CONDITIONS",
            "why_not_pass": "nothing was demonstrated",
            "why_not_fail": (
                "FAIL would assert the outputs were tested and found wanting. They were not "
                "tested. FAIL is false in the opposite direction, and is as misleading as PASS."
            ),
            "implementation_requirement": (
                "Compute the applicable set first, branch on len(applicable) == 0 before any "
                "aggregation, and emit the third state. Do not reach the aggregation at all."
            ),
        },

        "the_part_that_matters_most": {
            "finding": "This trap was CREATED by the fix for the earlier ones.",
            "before_the_fix": (
                "Without honest not_applicable, a gate in this state would have reported a "
                "misleading met/total - counting conditions that could not fail as if they had "
                "been satisfied."
            ),
            "after_the_fix": (
                "With honest not_applicable, the same gate can report a spurious PASS, because "
                "the denominator correctly empties and the aggregation vacuously succeeds."
            ),
            "the_lesson": (
                "A correctness fix moved the failure rather than removing it, and moved it "
                "somewhere QUIETER. The earlier bug was visible in a number a reader could "
                "question - 5 of 5 on a deterministic fit looks odd. The new bug is invisible: "
                "every individual condition behaves correctly and the defect exists only in the "
                "aggregation over them."
            ),
            "generalised": (
                "When a fix changes how results are aggregated or counted, re-audit the "
                "aggregation itself, not only the components. A component-level correctness fix "
                "can introduce a composite-level defect that no component-level test will catch."
            ),
            "related_prior_lesson": (
                "This is the same shape as the W2 finding that per-artifact checks cannot see a "
                "join: correctness at the part does not imply correctness at the whole."
            ),
        },

        "the_vacuity_shape_five_instances": [
            {"n": 1, "instance": "seed_direction_consistency on a deterministic fit", "scope": "condition", "caught": "at design time", "would_have_paid": "5-of-5 by construction"},
            {"n": 2, "instance": "qualifying_lineage_effect_floor scored while also filtering upstream", "scope": "condition", "caught": "before any result was observed", "would_have_paid": "always true for anything reaching it"},
            {"n": 3, "instance": "qualifying_set_dominance_guard on a closed composition", "scope": "condition", "caught": "only after the result, by the closure diagnostic", "would_have_paid": "satisfied by the dominant component's arithmetic complement", "note": "the one that got through"},
            {"n": 4, "instance": "a gate condition the evaluator never computed", "scope": "evaluator", "caught": "at review", "would_have_paid": "reports as absent rather than failing"},
            {"n": 5, "instance": "zero applicable conditions aggregated with all()", "scope": "GATE", "caught": "while drafting a gate for a lane that cannot run", "would_have_paid": "PASS having tested nothing", "note": "created by the fix for 1 and 2; the quietest of the five"},
        ],

        "do_not_retrofit": {
            "instruction": "Do NOT open a v3 of the W2 gate to retrofit this rule.",
            "w2_gate_path": str(a.w2_gate),
            "w2_gate_sha256": sha256_file(a.w2_gate),
            "reason": (
                "That gate is frozen forward-only and its completed run produced a real verdict "
                "with a non-empty applicable set, so the rule could not have changed its "
                "outcome. Churning it would cost more than it buys."
            ),
            "applies_to": "gates registered after this record",
        },

        "read_only_inputs_unmodified": {
            "w3_record_sha256": sha256_file(a.w3_record / "w3_recorded_negative.toml"),
            "typed_evidence_graph_sha256": sha256_file(a.typed_evidence_graph),
        },
        "external_development_only": True,
        "champion_eligible": False,
    }
    a.output.mkdir(parents=True)
    (a.output / "zero_applicable_gate_convention.json").write_text(
        json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(a.output, {
        "artifact_class": "campaign_wide_gate_convention",
        "convention": "zero_applicable_conditions",
        "scope": "campaign_wide",
        "modifies_nothing": True,
        "status": "passed",
    })
    print(json.dumps({"rule": rec["the_rule"]["required_third_state"],
                      "instances": len(rec["the_vacuity_shape_five_instances"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
