#!/usr/bin/env python3
"""Independently rederive the sealed blinded Stage-A QC decisions."""

from __future__ import annotations

import json
from collections import Counter

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def bh(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda index: (p_values[index], index))
    result = [1.0] * n
    running = 1.0
    for position in range(n - 1, -1, -1):
        index = order[position]
        running = min(running, p_values[index] * n / (position + 1))
        result[index] = min(1.0, running)
    return result


def metric_pass(value: float, operator: str, threshold: float) -> bool:
    if operator == "ge":
        return value >= threshold
    if operator == "le":
        return value <= threshold
    if operator == "eq":
        return value == threshold
    raise RuntimeError(f"Unsupported QC operator: {operator}")


def main() -> None:
    seal_path = CANDIDATE_ROOT / "BLINDED_STAGE_A_QC_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") not in {
        "blinded_stage_a_qc_sealed", "blinded_stage_a_qc_sealed_no_samples_pass"
    }:
        raise RuntimeError("Invalid blinded Stage-A QC status")
    for field in ["sample_exclusions_frozen", "guide_response_exclusions_frozen"]:
        if seal.get(field) is not True:
            raise RuntimeError(f"Blinded QC did not freeze {field}")
    for field in [
        "scientific_condition_labels_opened", "scientific_outcomes_inspected",
        "unblinding_authorized", "stage_b_design_frozen",
    ]:
        if seal.get(field) is not False:
            raise RuntimeError(f"Blinded QC improperly sets {field}")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Blinded-QC output hash mismatch: {name}")

    manifest = read_tsv(CANDIDATE_ROOT / "blinded_qc_input_manifest.tsv")
    expected_roles = {
        "blinded_sample_file_manifest", "blinded_sample_qc_metrics",
        "target_independent_guide_response",
        "blinded_qc_signoff", "stage_a_execution_seal",
        "qc_unblinding_contract_seal", "qc_threshold_seal",
    }
    if {row["role"] for row in manifest} != expected_roles:
        raise RuntimeError("Blinded-QC input-manifest universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Blinded-QC source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Blinded-QC source hash drift: {path}")

    by_role = {row["role"]: PROJECT_ROOT / row["source_path"] for row in manifest}
    contract_row = next(row for row in manifest if row["role"] == "qc_unblinding_contract_seal")
    contract_seal = json.loads((PROJECT_ROOT / contract_row["source_path"]).read_text(encoding="utf-8"))
    q_threshold = float(contract_seal["guide_response_q_threshold"])
    effect_threshold = float(contract_seal["guide_response_abs_effect_threshold"])

    status = read_tsv(CANDIDATE_ROOT / "blinded_sample_qc_status.tsv")
    if len(status) != int(seal["n_randomized_samples"]):
        raise RuntimeError("Blinded-QC sample count drift")
    if len({row["sample_id"] for row in status}) != len(status):
        raise RuntimeError("Blinded-QC sample IDs are duplicated")
    n_pass = sum(yes(row["technical_qc_pass"]) for row in status)
    if n_pass != int(seal["n_samples_included"]):
        raise RuntimeError("Blinded-QC included-sample count drift")
    execution_root = by_role["stage_a_execution_seal"].parent
    randomization = read_tsv(execution_root / "frozen_randomization_manifest.tsv")
    expected_sample = {
        row["sample_id"]: (row["blinded_label"], row["assay_id"])
        for row in randomization
    }
    if len(expected_sample) != len(randomization):
        raise RuntimeError("Bound execution freeze has duplicate sample IDs")
    threshold_root = by_role["qc_threshold_seal"].parent
    thresholds = read_tsv(threshold_root / "frozen_blinded_assay_qc_thresholds.tsv")
    threshold_by_key = {
        (row["assay_id"], row["metric_name"]): row for row in thresholds
    }
    if len(threshold_by_key) != len(thresholds):
        raise RuntimeError("Independent threshold universe contains duplicates")
    metrics = read_tsv(by_role["blinded_sample_qc_metrics"])
    metric_by_key = {(row["sample_id"], row["metric_name"]): row for row in metrics}
    if len(metric_by_key) != len(metrics):
        raise RuntimeError("Independent blinded metric universe contains duplicates")
    expected_metrics = {
        (sample_id, metric_name)
        for sample_id, (_, assay_id) in expected_sample.items()
        for threshold_assay, metric_name in threshold_by_key
        if threshold_assay == assay_id
    }
    if set(metric_by_key) != expected_metrics:
        raise RuntimeError("Independent blinded metric universe differs from frozen design")
    observed_status = {row["sample_id"]: row for row in status}
    if set(observed_status) != set(expected_sample):
        raise RuntimeError("Blinded-QC status differs from frozen randomized sample universe")
    for sample_id, (_, assay_id) in expected_sample.items():
        failures = []
        for key, metric in metric_by_key.items():
            if key[0] != sample_id:
                continue
            threshold = threshold_by_key[(assay_id, key[1])]
            if not metric_pass(
                float(metric["metric_value"]), threshold["operator"],
                float(threshold["threshold"]),
            ):
                failures.append(key[1])
        row = observed_status[sample_id]
        expected_pass = not failures
        expected_status = (
            "included_after_blinded_qc" if expected_pass
            else "excluded_before_unblinding"
        )
        expected_reason = (
            "" if expected_pass
            else "failed_frozen_metrics:" + ";".join(sorted(failures))
        )
        if (
            yes(row["technical_qc_pass"]) != expected_pass
            or row["preunblinding_exclusion_status"] != expected_status
            or row["preunblinding_exclusion_reason"] != expected_reason
        ):
            raise RuntimeError("Independent blinded-QC exclusion rederivation failed")
        if yes(row["scientific_condition_inspected"]):
            raise RuntimeError("Blinded-QC status reports condition access")

    guide = read_tsv(CANDIDATE_ROOT / "guide_response_exclusion_genes.tsv")
    raw_guide = read_tsv(by_role["target_independent_guide_response"])
    adjusted = bh([float(row["p_value"]) for row in raw_guide])
    expected_guide = {
        row["gene_id"]: (
            q_value,
            q_value < q_threshold and abs(float(row["effect"])) >= effect_threshold,
        )
        for row, q_value in zip(raw_guide, adjusted)
    }
    if {row["gene_id"] for row in guide} != set(expected_guide):
        raise RuntimeError("Guide-response output universe differs from raw frozen input")
    for row in guide:
        q_value, expected = expected_guide[row["gene_id"]]
        if abs(float(row["q_value"]) - q_value) > 1e-8:
            raise RuntimeError("Independent guide-response BH rederivation failed")
        if yes(row["exclude_from_all_stage_a_axes"]) != expected:
            raise RuntimeError("Independent guide-response exclusion rederivation failed")
    if len(guide) != int(seal["n_guide_response_genes_tested"]):
        raise RuntimeError("Guide-response tested-universe count drift")
    if sum(yes(row["exclude_from_all_stage_a_axes"]) for row in guide) != int(seal["n_guide_response_genes_excluded"]):
        raise RuntimeError("Guide-response exclusion count drift")

    retention = read_tsv(CANDIDATE_ROOT / "blinded_qc_design_retention.tsv")
    counts = Counter(row["assay_id"] for row in status)
    included = Counter(
        row["assay_id"] for row in status if yes(row["technical_qc_pass"])
    )
    if {
        row["assay_id"]: (
            int(row["n_randomized_samples"]), int(row["n_included_after_blinded_qc"])
        )
        for row in retention
    } != {assay: (counts[assay], included[assay]) for assay in counts}:
        raise RuntimeError("Blinded-QC assay retention does not rederive")
    print(
        "BLINDED_STAGE_A_QC_VALIDATION_PASS "
        f"samples={len(status)} included={n_pass} guide_genes={len(guide)} "
        "conditions_opened=false outcomes_opened=false"
    )


if __name__ == "__main__":
    main()
