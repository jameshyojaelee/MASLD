#!/usr/bin/env python3
"""Re-derive every targeted pair export and freeze orientation evidence.

This stage reports pass/fail for every prespecified locus.  It does not choose
an experimental target: lineage observability, element accessibility, and
editability remain separate downstream gates.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
from collections import defaultdict
from pathlib import Path

from relay_common import CANDIDATE_ROOT, atomic_write_json, read_tsv, sha256_file, write_tsv


def number(value: object) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"Non-finite orientation value: {value!r}")
    return result


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def pair_family_verdict(rows: list[dict[str, object]]) -> tuple[bool, str, str]:
    passing = [row for row in rows if row["orientation_gate_pass"] == "true"]
    directions = {str(row["oriented_risk_effect"]) for row in passing}
    all_pairs_pass = len(passing) == len(rows)
    concordant = all_pairs_pass and len(directions) == 1
    reason = (
        "all_eligible_pairs_oriented_and_concordant"
        if concordant
        else "one_or_more_eligible_pairs_failed_orientation"
        if not all_pairs_pass
        else "eligible_pairs_have_conflicting_risk_expression_directions"
    )
    direction = next(iter(directions)) if concordant else "unresolved"
    return concordant, reason, direction


def main() -> None:
    worklist_path = CANDIDATE_ROOT / "orientation_worklist.tsv"
    worklist_seal_path = CANDIDATE_ROOT / "ORIENTATION_WORKLIST_SEALED.json"
    release_path = CANDIDATE_ROOT / "ORIENTATION_RELEASE.json"
    if release_path.exists():
        raise RuntimeError(f"Refusing to overwrite orientation release: {release_path}")
    for path in [worklist_path, worklist_seal_path]:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing frozen worklist input: {path}")
    worklist_seal = json.loads(worklist_seal_path.read_text(encoding="utf-8"))
    if sha256_file(worklist_path) != worklist_seal["orientation_worklist_sha256"]:
        raise RuntimeError("Orientation worklist differs from its pre-export seal")

    worklist = read_tsv(worklist_path)
    orientations: list[dict[str, object]] = []
    manifest: list[dict[str, object]] = []
    for task in worklist:
        paths = {
            "pair_summary": CANDIDATE_ROOT / task["required_pair_summary"],
            "variant_posterior": CANDIDATE_ROOT / task["required_variant_posterior"],
            "allele_audit": CANDIDATE_ROOT / task["required_allele_audit"],
        }
        for role, path in paths.items():
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"Missing {role} for {task['orientation_uid']}: {path}")
            manifest.append(
                {
                    "orientation_uid": task["orientation_uid"],
                    "role": role,
                    "relative_path": str(path.relative_to(CANDIDATE_ROOT)),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )

        pair_rows = read_tsv(paths["pair_summary"])
        audit_rows = read_tsv(paths["allele_audit"])
        variants = read_tsv(paths["variant_posterior"])
        if len(pair_rows) != 1 or len(audit_rows) != 1 or not variants:
            raise RuntimeError(f"Invalid pair-export cardinality: {task['orientation_uid']}")
        pair = pair_rows[0]
        audit = audit_rows[0]
        if pair["orientation_uid"] != task["orientation_uid"] or audit[
            "orientation_uid"
        ] != task["orientation_uid"]:
            raise RuntimeError(f"Pair-export UID mismatch: {task['orientation_uid']}")
        if pair["target_freeze_status"] != "prohibited_until_lineage_and_editability_gates":
            raise RuntimeError(f"Pair export opened target gate: {task['orientation_uid']}")
        if yes(audit["target_list_frozen"]):
            raise RuntimeError(f"Allele audit claims frozen targets: {task['orientation_uid']}")

        pp4_delta = number(pair["pp4_delta"])
        if abs(pp4_delta) > 1e-8 or audit["pair_pp4_reproduction_status"] != "pass":
            raise RuntimeError(f"Pair PP.H4 reproduction failed: {task['orientation_uid']}")
        posterior = [number(row["shared_posterior"]) for row in variants]
        posterior_sum = sum(posterior)
        if abs(posterior_sum - 1.0) > 1e-6:
            raise RuntimeError(
                f"Variant shared posterior does not sum to one for {task['orientation_uid']}: "
                f"{posterior_sum}"
            )
        signs = [int(float(row["risk_to_expression_sign"])) for row in variants]
        if any(sign not in {-1, 0, 1} for sign in signs):
            raise RuntimeError(f"Invalid orientation sign: {task['orientation_uid']}")
        positive = sum(weight for weight, sign in zip(posterior, signs) if sign > 0)
        negative = sum(weight for weight, sign in zip(posterior, signs) if sign < 0)
        consensus = max(positive, negative)
        direction = (
            "risk_increases_expression"
            if positive > negative
            else "risk_decreases_expression"
            if negative > positive
            else "unresolved"
        )
        if abs(consensus - number(pair["orientation_consensus"])) > 1e-8:
            raise RuntimeError(f"Orientation consensus mismatch: {task['orientation_uid']}")
        if direction != pair["oriented_risk_effect"]:
            raise RuntimeError(f"Orientation direction mismatch: {task['orientation_uid']}")

        ordered = sorted(
            zip(variants, posterior, signs),
            key=lambda item: (-item[1], item[0]["snp"]),
        )
        leave = ordered[1:]
        leave_total = sum(item[1] for item in leave)
        if leave_total > 0:
            leave_positive = sum(item[1] for item in leave if item[2] > 0) / leave_total
            leave_negative = sum(item[1] for item in leave if item[2] < 0) / leave_total
            leave_consensus = max(leave_positive, leave_negative)
            leave_direction = (
                "risk_increases_expression"
                if leave_positive > leave_negative
                else "risk_decreases_expression"
                if leave_negative > leave_positive
                else "unresolved"
            )
        else:
            leave_consensus = float("nan")
            leave_direction = "unresolved"
        recorded_leave = pair["leave_top_orientation_consensus"].strip()
        if recorded_leave:
            if not math.isfinite(leave_consensus) or abs(leave_consensus - number(recorded_leave)) > 1e-8:
                raise RuntimeError(f"Leave-top consensus mismatch: {task['orientation_uid']}")
        elif math.isfinite(leave_consensus):
            raise RuntimeError(f"Missing finite leave-top consensus: {task['orientation_uid']}")
        if leave_direction != pair["leave_top_oriented_risk_effect"]:
            raise RuntimeError(f"Leave-top direction mismatch: {task['orientation_uid']}")

        shared_95 = [row for row in variants if yes(row["in_shared_95"])]
        if not shared_95:
            raise RuntimeError(f"Empty shared 95% set: {task['orientation_uid']}")
        allele_resolved = all(
            row["gwas_effect_allele"]
            and row["gwas_other_allele"]
            and row["eqtl_effect_allele_source"]
            and row["eqtl_other_allele_source"]
            and (
                not yes(row["is_palindromic"])
                or row["gwas_strand"] == "same"
            )
            for row in shared_95
        )
        orientation_pass = (
            number(task["susie_pp4"]) > 0.5
            and consensus >= 0.80
            and leave_direction == direction
            and direction != "unresolved"
            and allele_resolved
        )
        orientations.append(
            {
                "orientation_uid": task["orientation_uid"],
                "pair_family_uid": task["pair_family_uid"],
                "coarse_locus_uid": task["coarse_locus_uid"],
                "gene_symbol": task["gene_symbol"],
                "ensembl_id": task["ensembl_id"],
                "gwas_name": task["gwas_name"],
                "trait": task["trait"],
                "tier": task["tier"],
                "phenotype_stratum": task["phenotype_stratum"],
                "source_dependence": task["source_dependence"],
                "gwas_evidence_family": task["gwas_evidence_family"],
                "gwas_independence_class": task["gwas_independence_class"],
                "counts_for_replication_breadth": task[
                    "counts_for_replication_breadth"
                ],
                "n_eligible_gwas_pairs_for_locus_gene": task[
                    "n_eligible_gwas_pairs_for_locus_gene"
                ],
                "susie_pp4": task["susie_pp4"],
                "selected_gwas_signal_index": pair["gwas_signal_index"],
                "selected_eqtl_signal_index": pair["eqtl_signal_index"],
                "selected_gwas_hit": pair["gwas_signal_hit"],
                "selected_eqtl_hit": pair["eqtl_signal_hit"],
                "orientation_consensus": consensus,
                "oriented_risk_effect": direction,
                "leave_top_orientation_consensus": (
                    "" if not math.isfinite(leave_consensus) else leave_consensus
                ),
                "leave_top_oriented_risk_effect": leave_direction,
                "credible_set_alleles_resolved": str(allele_resolved).lower(),
                "orientation_gate_pass": str(orientation_pass).lower(),
                "orientation_gate_reason": (
                    "pass_pending_lineage_and_editability"
                    if orientation_pass
                    else "failed_consensus_leave_top_or_allele_gate"
                ),
                "target_freeze_status": "prohibited_until_lineage_and_editability_gates",
            }
        )

    by_family: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in orientations:
        by_family[str(row["pair_family_uid"])].append(row)
    family_rows: list[dict[str, object]] = []
    for family_uid, rows in sorted(by_family.items()):
        passing = [row for row in rows if row["orientation_gate_pass"] == "true"]
        concordant, reason, family_direction = pair_family_verdict(rows)
        for row in rows:
            row["pair_family_orientation_concordant"] = str(concordant).lower()
            row["pair_family_gate_reason"] = reason
        family_rows.append(
            {
                "pair_family_uid": family_uid,
                "coarse_locus_uid": rows[0]["coarse_locus_uid"],
                "gene_symbol": rows[0]["gene_symbol"],
                "ensembl_id": rows[0]["ensembl_id"],
                "n_eligible_pairs": len(rows),
                "n_pair_orientation_pass": len(passing),
                "n_distinct_gwas": len({str(row["gwas_name"]) for row in rows}),
                "n_distinct_gwas_evidence_families": len(
                    {str(row["gwas_evidence_family"]) for row in rows}
                ),
                "n_breadth_eligible_gwas_evidence_families": len(
                    {
                        str(row["gwas_evidence_family"])
                        for row in rows
                        if yes(row["counts_for_replication_breadth"])
                    }
                ),
                "n_distinct_phenotype_strata": len(
                    {str(row["phenotype_stratum"]) for row in rows}
                ),
                "family_oriented_risk_effect": (
                    family_direction
                ),
                "pair_family_orientation_concordant": str(concordant).lower(),
                "pair_family_gate_reason": reason,
                "independence_interpretation": (
                    "multiple_GWAS_are_directional_support_not_automatically_"
                    "independent_replication"
                ),
                "target_freeze_status": (
                    "prohibited_until_lineage_accessibility_and_guideability"
                ),
            }
        )

    orientation_path = CANDIDATE_ROOT / "shared_signal_orientation.tsv"
    family_path = CANDIDATE_ROOT / "shared_signal_orientation_families.tsv"
    manifest_path = CANDIDATE_ROOT / "pair_export_manifest.tsv"
    orientation_fields = [
        "orientation_uid",
        "pair_family_uid",
        "coarse_locus_uid",
        "gene_symbol",
        "ensembl_id",
        "gwas_name",
        "trait",
        "tier",
        "phenotype_stratum",
        "source_dependence",
        "gwas_evidence_family",
        "gwas_independence_class",
        "counts_for_replication_breadth",
        "n_eligible_gwas_pairs_for_locus_gene",
        "susie_pp4",
        "selected_gwas_signal_index",
        "selected_eqtl_signal_index",
        "selected_gwas_hit",
        "selected_eqtl_hit",
        "orientation_consensus",
        "oriented_risk_effect",
        "leave_top_orientation_consensus",
        "leave_top_oriented_risk_effect",
        "credible_set_alleles_resolved",
        "orientation_gate_pass",
        "orientation_gate_reason",
        "pair_family_orientation_concordant",
        "pair_family_gate_reason",
        "target_freeze_status",
    ]
    write_tsv(orientation_path, orientations, orientation_fields)
    write_tsv(family_path, family_rows, list(family_rows[0]))
    write_tsv(
        manifest_path,
        manifest,
        ["orientation_uid", "role", "relative_path", "size_bytes", "sha256"],
    )
    payload = {
        "status": "orientation_evidence_frozen_targets_not_selected",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "worklist_seal_sha256": sha256_file(worklist_seal_path),
        "worklist_sha256": sha256_file(worklist_path),
        "n_prespecified_pairs": len(worklist),
        "n_pair_families": len(family_rows),
        "n_orientation_pair_gate_pass": sum(
            row["orientation_gate_pass"] == "true" for row in orientations
        ),
        "n_pair_family_orientation_concordant": sum(
            row["pair_family_orientation_concordant"] == "true"
            for row in family_rows
        ),
        "scientific_outcomes_inspected": False,
        "experimental_targets_frozen": False,
        "next_gate": "source-lineage expression, element accessibility, guideability, and exact-edit feasibility",
        "shared_signal_orientation_sha256": sha256_file(orientation_path),
        "shared_signal_orientation_families_sha256": sha256_file(family_path),
        "pair_export_manifest_sha256": sha256_file(manifest_path),
    }
    atomic_write_json(release_path, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
