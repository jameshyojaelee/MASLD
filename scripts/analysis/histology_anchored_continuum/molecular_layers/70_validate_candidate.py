#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_finite_number(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def command_output(command: list[str]) -> str:
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    return result.stdout


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: 70_validate_candidate.py CANDIDATE_ROOT")
    root = Path(sys.argv[1]).resolve()
    if not root.is_dir():
        raise SystemExit(f"candidate root is absent: {root}")

    contract_path = root / "inputs" / "frozen_contract.json"
    contract = json.loads(contract_path.read_text())
    checks: list[dict[str, object]] = []

    def check(check_id: str, passed: bool, observed: object, expected: object) -> None:
        checks.append(
            {
                "check_id": check_id,
                "pass": bool(passed),
                "observed": observed,
                "expected": expected,
            }
        )

    required = [
        "inputs/input_manifest.tsv",
        "inputs/analysis_code_manifest.tsv",
        "inputs/source_overlap_matrix.tsv",
        "bulk/continuum_gene_cohort.tsv.gz",
        "bulk/continuum_gene_meta.tsv.gz",
        "bulk/continuum_membership.tsv.gz",
        "bulk/signature_loo_gene_meta.tsv",
        "pathway_tf/pathways/collection_index.tsv",
        "pathway_tf/tf/continuum_tf_activity.tsv",
        "pathway_tf/tf/regulon_meta_analysis.tsv",
        "programs/nmf/nmf_continuum_membership.tsv",
        "programs/hotspot/hotspot_continuum_membership.tsv",
        "programs/hotspot/hotspot_all117_within_stage_permutations.tsv",
        "paired/paired_participant_manifest.tsv",
        "paired/paired_program_results.tsv",
        "decision/continuum_membership_registry.tsv.gz",
        "decision/integration_gates.tsv",
        "decision/candidate_decision.json",
        "figures/panels/fig4f_continuum_program_trajectories.pdf",
        "figures/figure_manifest.tsv",
        "figures/review_assembly_manifest.tsv",
    ]
    missing = [relative for relative in required if not (root / relative).is_file()]
    check("required_artifacts", not missing, ";".join(missing), "none missing")

    input_manifest = read_tsv(root / "inputs" / "input_manifest.tsv")
    input_hashes_valid = all(
        Path(row["path"]).is_file() and sha256(Path(row["path"])) == row["sha256"]
        for row in input_manifest
    )
    check("input_hashes_recompute", input_hashes_valid, input_hashes_valid, True)
    code_manifest = read_tsv(root / "inputs" / "analysis_code_manifest.tsv")
    code_hashes_valid = all(
        Path(row["path"]).is_file() and sha256(Path(row["path"])) == row["sha256"]
        for row in code_manifest
    )
    check("analysis_code_hashes_recompute", code_hashes_valid,
          code_hashes_valid, True)

    signature = read_tsv(root / "bulk" / "signature_gene_ledger.tsv")
    observed_signature = sum(row.get("observed_in_F_five", "").upper() == "TRUE" for row in signature)
    check("signature_published_145", len(signature) == 145, len(signature), 145)
    check("signature_observed_139", observed_signature == 139, observed_signature, 139)

    transcript = read_tsv(root / "bulk" / "continuum_membership.tsv.gz")
    check("transcript_family_23370", len(transcript) == 23370, len(transcript), 23370)
    axis_members = sum(row.get("axis_component", "").upper() == "TRUE" for row in transcript)
    check("signature_axes_not_ordinary_members", axis_members == 139, axis_members, 139)

    loo = read_tsv(root / "bulk" / "signature_loo_gene_meta.tsv")
    check("loo_family_139_by_2", len(loo) == 278, len(loo), 278)

    for collection, family_size in contract["pathway_collections"].items():
        meta = read_tsv(root / "pathway_tf" / "pathways" / collection / "meta_analysis.tsv")
        ids = {row["set_id"] for row in meta}
        check(
            f"pathway_{collection}_complete",
            len(ids) == family_size and len(meta) == family_size * 2,
            f"{len(ids)} ids/{len(meta)} rows",
            f"{family_size} ids/{family_size * 2} rows",
        )

    tf = read_tsv(root / "pathway_tf" / "tf" / "continuum_tf_activity.tsv")
    tf_ids = {row["tf"] for row in tf}
    check("tf_family_268", len(tf_ids) == 268 and len(tf) == 268 * 4,
          f"{len(tf_ids)} ids/{len(tf)} rows", "268 ids/1072 rows")
    regulon_testability = read_tsv(root / "pathway_tf" / "tf" / "regulon_testability.tsv")
    display_tfs = set(contract["display_tfs"])
    evaluation_cohorts = set(contract["evaluation_cohorts"])
    display_status = [
        row for row in regulon_testability
        if row["tf"] in display_tfs and row["dataset"] in evaluation_cohorts
    ]
    check("display_regulon_testability_explicit", len(display_status) == 8,
          len(display_status), 8)
    testable_display = {
        (row["tf"], row["dataset"])
        for row in display_status if row["testable"].upper() == "TRUE"
    }
    regulon_scores = read_tsv(root / "pathway_tf" / "tf" / "donor_regulon_scores.tsv.gz")
    finite_display_scores = {
        (row["tf"], row["dataset"])
        for row in regulon_scores
        if row["tf"] in display_tfs and row["dataset"] in evaluation_cohorts
        and is_finite_number(row["regulon_score"])
    }
    missing_finite_regulons = sorted(testable_display - finite_display_scores)
    check("testable_display_regulons_have_finite_scores", not missing_finite_regulons,
          ";".join(f"{tf_name}:{cohort}" for tf_name, cohort in missing_finite_regulons),
          "none")

    nmf = read_tsv(root / "programs" / "nmf" / "nmf_continuum_membership.tsv")
    check("nmf_family_10", len(nmf) == 10, len(nmf), 10)

    hotspot = read_tsv(root / "programs" / "hotspot" / "hotspot_continuum_membership.tsv")
    check("hotspot_family_117", len(hotspot) == 117, len(hotspot), 117)
    hotspot_testability = read_tsv(root / "programs" / "hotspot" / "hotspot_testability.tsv")
    pdgfra = [row for row in hotspot_testability if "PDGFRA" in row.get("module_name", "")]
    pdgfra_untestable = len(pdgfra) == 2 and all(
        row.get("testable", "").upper() == "FALSE" for row in pdgfra
    )
    check("pdgfra_program_explicit_untestable", pdgfra_untestable, len(pdgfra), 2)

    permutations = read_tsv(
        root / "programs" / "hotspot" / "hotspot_all117_within_stage_permutations.tsv"
    )
    replicates = {row.get("permutation_replicates") for row in permutations}
    check("hotspot_permutations_10000", replicates == {"10000"}, sorted(replicates), ["10000"])
    check("hotspot_permutation_family_rows", len(permutations) == 117 * 4,
          len(permutations), 468)

    paired = read_tsv(root / "paired" / "paired_participant_manifest.tsv")
    donor_ids = {row["donor_id"] for row in paired}
    check("paired_54_participants_108_biopsies", len(donor_ids) == 54 and len(paired) == 108,
          f"{len(donor_ids)} donors/{len(paired)} rows", "54 donors/108 rows")

    gates = read_tsv(root / "decision" / "integration_gates.tsv")
    blocking_failures = [
        row["check_id"] for row in gates
        if row.get("blocking", "").upper() == "TRUE" and row.get("pass", "").upper() != "TRUE"
    ]
    check("integration_blocking_gates", not blocking_failures,
          ";".join(blocking_failures), "none")

    decision = json.loads((root / "decision" / "candidate_decision.json").read_text())
    check("candidate_not_auto_promoted", decision.get("automatic_promotion") is False,
          decision.get("automatic_promotion"), False)
    check("figure3_unchanged", decision.get("current_figure3_modified") is False,
          decision.get("current_figure3_modified"), False)

    pdfs = sorted((root / "figures").rglob("*.pdf"))
    check("candidate_pdfs_exist", len(pdfs) >= 16, len(pdfs), ">=16")
    pdf_failures: list[str] = []
    for pdf in pdfs:
        info = command_output(["pdfinfo", str(pdf)])
        pages = next((line.split(":", 1)[1].strip() for line in info.splitlines()
                      if line.startswith("Pages:")), "")
        if pages != "1":
            pdf_failures.append(f"{pdf.name}:pages={pages}")
        fonts = command_output([
            "gs", "-q", "-dNOSAFER", "-dPDFDEBUG", "-sDEVICE=nullpage",
            "-o", "/dev/null", str(pdf),
        ])
        if "/Subtype /Type3" in fonts or "Dingbats" in fonts:
            pdf_failures.append(f"{pdf.name}:forbidden_font")
        if "Helvetica" not in fonts:
            pdf_failures.append(f"{pdf.name}:helvetica_missing")
        if "Helvetica-Bold" in fonts:
            pdf_failures.append(f"{pdf.name}:bold_font")
    check("pdf_one_page_helvetica_no_bold_type3_or_dingbats", not pdf_failures,
          ";".join(pdf_failures), "none")

    figure_manifest = read_tsv(root / "figures" / "figure_manifest.tsv")
    manifest_paths_exist = all(Path(row["path"]).is_file() for row in figure_manifest)
    check("figure_manifest_paths_exist", manifest_paths_exist,
          manifest_paths_exist, True)
    source_tables = list((root / "figures" / "source_tables").glob("*.tsv*"))
    check("figure_source_tables_present", len(source_tables) >= 16,
          len(source_tables), ">=16")

    forbidden_figure3 = [str(path.relative_to(root)) for path in root.rglob("*")
                         if path.is_file() and "figure3" in path.name.lower()]
    check("no_figure3_outputs", not forbidden_figure3,
          ";".join(forbidden_figure3), "none")

    validation_dir = root / "validation"
    validation_dir.mkdir(parents=True, exist_ok=False)
    manifest_paths = sorted(path for path in root.rglob("*") if path.is_file())
    with (validation_dir / "candidate_file_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        for path in manifest_paths:
            writer.writerow([path.relative_to(root), path.stat().st_size, sha256(path)])

    with (validation_dir / "validation_checks.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["check_id", "pass", "observed", "expected"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(checks)

    passed = all(bool(row["pass"]) for row in checks)
    summary = {
        "candidate_root": str(root),
        "validation_pass": passed,
        "n_checks": len(checks),
        "n_failed": sum(not bool(row["pass"]) for row in checks),
        "failed_checks": [row["check_id"] for row in checks if not bool(row["pass"])],
        "release_state": "candidate_only",
    }
    (validation_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
