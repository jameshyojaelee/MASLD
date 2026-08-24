#!/usr/bin/env python3
"""Independent validation and readiness sealing for promoted-COLOC context."""

from __future__ import annotations

import csv
import hashlib
import math
import os
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import pysam


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = (ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID).resolve()


def read(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def classify(mapped: float, shared: float, left: float, right: float) -> str:
    if mapped < 0.95:
        return "untestable"
    if shared >= 0.5:
        return "replicated_accessible"
    if max(left, right) >= 0.5:
        return "source_dependent"
    if max(shared, left, right) > 0:
        return "partial"
    return "indeterminate"


def atomic(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    if path.exists() or path.is_symlink():
        raise RuntimeError(f"refusing to overwrite genetics validator output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, encoding="utf-8", newline="") as handle:
        temp = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def main() -> None:
    if not (CANDIDATE / "NON_GENETIC_READY").is_file():
        raise RuntimeError("NON_GENETIC_READY is required before genetics sealing")
    if not (CANDIDATE / "CANDIDATE_AUDIT_READY").is_file():
        raise RuntimeError("CANDIDATE_AUDIT_READY is required before genetics sealing")
    prepared = CANDIDATE / "genetics/prepared"
    if not (prepared / "PROMOTED_UPSTREAM_GATE.tsv").is_file():
        raise RuntimeError("promoted COLOC gate is absent")
    out = CANDIDATE / "genetics/validation"
    if (
        out.exists() or out.is_symlink()
        or (CANDIDATE / "GENETIC_READY").exists()
        or (CANDIDATE / "GENETIC_READY").is_symlink()
        or (CANDIDATE / "FULL_READY").exists()
        or (CANDIDATE / "FULL_READY").is_symlink()
    ):
        raise RuntimeError("refusing to overwrite genetics readiness outputs")
    gate_rows = read(prepared / "PROMOTED_UPSTREAM_GATE.tsv")
    if len(gate_rows) != 1 or gate_rows[0]["status"] != "PROMOTED":
        raise RuntimeError("promoted upstream gate is invalid")
    replay_inputs = read(prepared / "replay_input_manifest.tsv")
    if not replay_inputs:
        raise RuntimeError("replay input manifest is empty")
    for row in replay_inputs:
        path = Path(row["path"])
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            raise RuntimeError(f"replay input hash mismatch: {path}")
    if sha256(prepared / "replay_input_manifest.tsv") != gate_rows[0]["replay_input_manifest_sha256"]:
        raise RuntimeError("replay input manifest hash disagrees with promoted gate")
    if sha256(prepared / "replay_plan.tsv") != gate_rows[0]["replay_plan_sha256"]:
        raise RuntimeError("replay plan hash disagrees with promoted gate")
    if sha256(prepared / "replay_runner/06_susie_coloc_replay.R") != gate_rows[0]["replay_runner_sha256"]:
        raise RuntimeError("isolated replay runner hash disagrees with promoted gate")
    batch_manifests = [
        CANDIDATE / f"genetics/replay_execution/batch_{batch}/batch_manifest.tsv"
        for batch in range(1, 6)
    ]
    if not all(path.is_file() and not path.is_symlink() for path in batch_manifests):
        raise RuntimeError("exactly five validated replay batch manifests are required")
    replay_artifacts: list[Path] = []
    for batch, manifest_path in enumerate(batch_manifests, start=1):
        batch_root = manifest_path.parent.resolve()
        rows = read(manifest_path)
        expected_batch_pairs = {
            (row["gwas_name"], row["ensembl"])
            for row in read(prepared / f"replay_batches/batch_{batch}.tsv")
        }
        observed_pair_artifacts: Counter[tuple[str, str]] = Counter()
        for row in rows:
            if row["release_id"] != RELEASE_ID or int(row["batch_id"]) != batch:
                raise RuntimeError("replay batch manifest identity mismatch")
            artifact = (batch_root / row["artifact"]).resolve()
            try:
                artifact.relative_to(batch_root)
            except ValueError as error:
                raise RuntimeError("replay artifact escapes its batch root") from error
            if (
                not artifact.is_file()
                or artifact.stat().st_size != int(row["bytes"])
                or sha256(artifact) != row["sha256"]
            ):
                raise RuntimeError(f"validated replay artifact changed: {artifact}")
            observed_pair_artifacts[(row["gwas_name"], row["ensembl"])] += 1
            replay_artifacts.append(artifact)
        if set(observed_pair_artifacts) != expected_batch_pairs or any(
            count != 2 for count in observed_pair_artifacts.values()
        ):
            raise RuntimeError(f"replay batch {batch} artifact family is incomplete")
    all_rows = read(CANDIDATE / "genetics/context/genetic_lineage_context_all_pairs.tsv")
    primary = read(CANDIDATE / "genetics/context/genetic_lineage_context_primary_pairs.tsv")
    if not all_rows or not primary:
        raise RuntimeError("genetic lineage context is empty")
    keys = Counter((row["gwas_name"], row["ensembl"], row["signal_pair_index"]) for row in all_rows)
    if any(value != 5 for value in keys.values()):
        raise RuntimeError("each signal pair must have exactly five lineage rows")
    primary_keys = Counter((row["gwas_name"], row["ensembl"]) for row in primary)
    if any(value != 5 for value in primary_keys.values()):
        raise RuntimeError("each gene-study pair must have one primary pair across five lineages")
    plan = read(prepared / "replay_plan.tsv")
    expected_pairs = {(row["gwas_name"], row["ensembl"]): float(row["promoted_pp_h4_susie"]) for row in plan}
    if set(primary_keys) != set(expected_pairs):
        raise RuntimeError("replay did not return every promoted supported gene-study pair")
    replay_max: dict[tuple[str, str], float] = {}
    for row in all_rows:
        key = (row["gwas_name"], row["ensembl"])
        replay_max[key] = max(replay_max.get(key, float("-inf")), float(row["pp_h4"]))
    for key, expected in expected_pairs.items():
        # Amended 2026-08-20 from rel_tol=1e-7 to 1e-4, completing amendment 002.
        # This is the second copy of the promoted-reproduction contract; the first
        # is 22_validate_replay_batch.py:130. Same justification and same evidence:
        # across all 816 validated pairs the median drift is 5.40e-12, max relative
        # 6.60e-05, max absolute 4.54e-05, and 0 pairs exceed 1e-4. The promoted pair
        # closest to the 0.5 selection cut is 0.503716, 82x the largest absolute
        # drift, so no gene-level membership call can move.
        # See REPLAY_TOLERANCE_AMENDMENT.md in the candidate root.
        if not math.isclose(replay_max[key], expected, rel_tol=1e-4, abs_tol=1e-9):
            raise RuntimeError(f"replayed maximum PP.H4 does not reproduce promoted aggregate: {key}")
    signal_metadata: dict[tuple[str, str], set[tuple[int, float, str, str]]] = {}
    for row in all_rows:
        key = (row["gwas_name"], row["ensembl"])
        signal_metadata.setdefault(key, set()).add((
            int(row["signal_pair_index"]), float(row["pp_h4"]),
            row["gwas_signal"], row["eqtl_signal"],
        ))
    selected_primary = {
        key: sorted(values, key=lambda value: (-value[1], value[2], value[3], value[0]))[0][0]
        for key, values in signal_metadata.items()
    }
    observed_primary = {
        key: {int(row["signal_pair_index"]) for row in primary
              if (row["gwas_name"], row["ensembl"]) == key}
        for key in expected_pairs
    }
    if any(indexes != {selected_primary[key]} for key, indexes in observed_primary.items()):
        raise RuntimeError("primary signal-pair selection is not deterministic")
    for row in all_rows:
        mapped = float(row["mapped_posterior_mass"])
        lost = float(row["lost_posterior_mass"])
        left = float(row["gse244832_accessible_mass"])
        right = float(row["gse281367_accessible_mass"])
        shared = float(row["shared_accessible_mass"])
        any_mass = float(row["any_accessible_mass"])
        pp4 = float(row["pp_h4"])
        if not math.isclose(mapped + lost, 1.0, abs_tol=1e-8):
            raise RuntimeError("mapped and lost posterior mass do not sum to one")
        if shared > min(left, right) + 1e-10 or any_mass + 1e-10 < max(left, right):
            raise RuntimeError("lineage posterior mass set relation failed")
        if not math.isclose(any_mass, left + right - shared, abs_tol=1e-10):
            raise RuntimeError("lineage posterior union mass mismatch")
        if not math.isclose(float(row["gse244832_unique_accessible_mass"]), left - shared, abs_tol=1e-10):
            raise RuntimeError("GSE244832 unique posterior mass mismatch")
        if not math.isclose(float(row["gse281367_unique_accessible_mass"]), right - shared, abs_tol=1e-10):
            raise RuntimeError("GSE281367 unique posterior mass mismatch")
        if not math.isclose(float(row["joint_any_accessible_mass"]), pp4 * any_mass, abs_tol=1e-10):
            raise RuntimeError("joint accessible mass mismatch")
        if not math.isclose(float(row["joint_shared_accessible_mass"]), pp4 * shared, abs_tol=1e-10):
            raise RuntimeError("joint shared accessible mass mismatch")
        if not math.isclose(float(row["joint_gse244832_accessible_mass"]), pp4 * left, abs_tol=1e-10):
            raise RuntimeError("joint GSE244832 accessible mass mismatch")
        if not math.isclose(float(row["joint_gse281367_accessible_mass"]), pp4 * right, abs_tol=1e-10):
            raise RuntimeError("joint GSE281367 accessible mass mismatch")
        if row["evidence_state"] != classify(mapped, shared, left, right):
            raise RuntimeError("genetic evidence state mismatch")
    audit = read(CANDIDATE / "genetics/context/variant_liftover_audit.tsv")
    liftover_inputs = read(CANDIDATE / "genetics/liftover/liftover_input_manifest.tsv")
    expected_liftover_roles = {"hg19_to_hg38_chain"} | {
        f"replay_batch_manifest_{batch}" for batch in range(1, 6)
    }
    if {row["role"] for row in liftover_inputs} != expected_liftover_roles:
        raise RuntimeError("liftover input manifest does not contain the complete input family")
    for row in liftover_inputs:
        path = Path(row["path"])
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise RuntimeError(f"liftover input changed: {path}")
    posterior_sums: Counter[tuple[str, str, str]] = Counter()
    for row in audit:
        posterior_sums[(row["gwas_name"], row["ensembl"], row["signal_pair_index"])] += float(row["snp_pp_h4"])
    if any(not math.isclose(value, 1.0, abs_tol=1e-6) for value in posterior_sums.values()):
        raise RuntimeError("pre-liftover signal-pair posterior does not sum to one")
    mapped_sums: Counter[tuple[str, str, str]] = Counter()
    frozen = read(CANDIDATE / "input_manifest.tsv")
    fasta_rows = [row for row in frozen if row["role"] == "hg38_fasta"]
    if len(fasta_rows) != 1:
        raise RuntimeError("frozen hg38 FASTA is not unique")
    fasta_path = Path(fasta_rows[0]["relative_or_absolute_path"])
    if not fasta_path.is_absolute():
        fasta_path = ROOT / fasta_path
    if sha256(fasta_path) != fasta_rows[0]["sha256"]:
        raise RuntimeError("frozen hg38 FASTA hash mismatch")
    fasta = pysam.FastaFile(str(fasta_path))
    for row in audit:
        key = (row["gwas_name"], row["ensembl"], row["signal_pair_index"])
        if row["mapping_status"] == "mapped":
            mapped_sums[key] += float(row["snp_pp_h4"])
            start0 = int(row["hg38_position_1based"]) - 1
            observed_ref = fasta.fetch(
                row["hg38_chrom"], start0, start0 + len(row["hg38_ref"])
            ).upper()
            if observed_ref != row["hg38_ref"].upper():
                raise RuntimeError("independent hg38 reference-allele validation failed")
    fasta.close()
    for key, total in posterior_sums.items():
        context_mapped = {
            float(row["mapped_posterior_mass"])
            for row in all_rows
            if (row["gwas_name"], row["ensembl"], row["signal_pair_index"]) == key
        }
        if len(context_mapped) != 1 or not math.isclose(
            next(iter(context_mapped)), mapped_sums[key], abs_tol=1e-8
        ):
            raise RuntimeError("mapped posterior mass does not reproduce liftover audit")
    figure_manifest = read(CANDIDATE / "genetics/figures/figure_manifest.tsv")
    if len(figure_manifest) != 1:
        raise RuntimeError("genetic Figure 4D candidate manifest must contain one panel")
    figure = figure_manifest[0]
    pdf = CANDIDATE / figure["pdf"]
    source = CANDIDATE / figure["source_table"]
    if sha256(pdf) != figure["pdf_sha256"] or sha256(source) != figure["source_sha256"]:
        raise RuntimeError("genetic Figure 4D checksum mismatch")
    pdfinfo = subprocess.run(["pdfinfo", str(pdf)], text=True, capture_output=True, check=True).stdout
    if "Pages:           1" not in pdfinfo:
        raise RuntimeError("genetic Figure 4D candidate is not a one-page PDF")
    size_match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdfinfo)
    # Page-size tolerance amended 2026-08-20 from 0.2 pt to 1.0 pt (amendment 004).
    # 13_render_genetic.R:110-111 renders through cairo_pdf whenever cairo is
    # available, and cairo truncates the PDF page box to whole points: a 5.4 x 2.8
    # inch request writes 388 x 201 pts, not 388.8 x 201.6. The gaps of 0.8 and 0.6
    # pt exceeded the old tolerance deterministically, so this check could never
    # pass on a cairo host. Switching the renderer to base pdf() would emit exact
    # fractional points but is the device prone to Type 3 fonts, which the check
    # immediately below rejects. One point is 1/72 inch, so 1.0 still catches a
    # genuinely wrong page size. The sealed non-genetic validator asserts no
    # numeric page dimensions at all (10_validate_non_genetic.py:375).
    if size_match is None or not (
        math.isclose(float(size_match.group(1)), 388.8, abs_tol=1.0)
        and math.isclose(float(size_match.group(2)), 201.6, abs_tol=1.0)
    ):
        raise RuntimeError("genetic Figure 4D candidate dimensions are not 5.4 x 2.8 inches")
    if b"/Subtype /Type3" in pdf.read_bytes():
        raise RuntimeError("genetic Figure 4D candidate contains a Type 3 font")
    figure_source = read(source)
    if len(figure_source) != 50:
        raise RuntimeError("Figure 4D source must complete 2 trait classes x 5 lineages x 5 states")
    expected_traits = {"direct_MASLD", "liver_enzyme"}
    expected_lineages = {"hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"}
    expected_states = {
        "replicated_accessible", "source_dependent", "partial", "indeterminate", "untestable"
    }
    source_cells = Counter(
        (row["trait_class"], row["lineage"], row["evidence_state"])
        for row in figure_source
    )
    if (
        {key[0] for key in source_cells} != expected_traits
        or {key[1] for key in source_cells} != expected_lineages
        or {key[2] for key in source_cells} != expected_states
        or any(count != 1 for count in source_cells.values())
    ):
        raise RuntimeError("Figure 4D source is not the complete prespecified state grid")
    source_totals: dict[tuple[str, str], int] = Counter()
    source_denominators: dict[tuple[str, str], set[int]] = {}
    for row in figure_source:
        key = (row["trait_class"], row["lineage"])
        source_totals[key] += int(row["n_pairs"])
        source_denominators.setdefault(key, set()).add(int(row["denominator"]))
    if len(source_totals) != 10 or any(
        len(source_denominators[key]) != 1
        or source_totals[key] != next(iter(source_denominators[key]))
        for key in source_totals
    ):
        raise RuntimeError("Figure 4D denominators do not reproduce state counts")
    primary_denominators = Counter(
        (row["trait_class"], row["gwas_name"], row["ensembl"])
        for row in primary if row["lineage"] == "hepatocyte"
    )
    expected_denominators = Counter(key[0] for key in primary_denominators)
    for trait in expected_traits:
        for lineage in expected_lineages:
            denominator = next(iter(source_denominators[(trait, lineage)]))
            if denominator != expected_denominators[trait]:
                raise RuntimeError("Figure 4D denominator does not match primary pairs")
    validation = [{
        "release_id": RELEASE_ID,
        "check": "variant_consistent_lineage_context",
        "status": "PASS",
        "n_signal_pairs": len(keys),
        "n_primary_gene_study_pairs": len(primary_keys),
        "n_variant_rows": len(audit),
    }]
    atomic(
        out / "validation_results.tsv",
        ("release_id", "check", "status", "n_signal_pairs", "n_primary_gene_study_pairs", "n_variant_rows"),
        validation,
    )
    artifacts = [
        prepared / "PROMOTED_UPSTREAM_GATE.tsv",
        prepared / "replay_input_manifest.tsv",
        prepared / "replay_plan.tsv",
        prepared / "replay_runner/06_susie_coloc_replay.R",
        prepared / "replay_runner/genetics_export_helpers.R",
        *[prepared / f"replay_batches/batch_{batch}.tsv" for batch in range(1, 6)],
        *batch_manifests,
        *replay_artifacts,
        CANDIDATE / "genetics/liftover/variant_liftover_raw.tsv.gz",
        CANDIDATE / "genetics/liftover/liftover_input_manifest.tsv",
        CANDIDATE / "genetics/liftover/sessionInfo.txt",
        CANDIDATE / "genetics/context/variant_liftover_audit.tsv",
        CANDIDATE / "genetics/context/genetic_lineage_context_all_pairs.tsv",
        CANDIDATE / "genetics/context/genetic_lineage_context_primary_pairs.tsv",
        CANDIDATE / "genetics/context/genetics_context_input_manifest.tsv",
        CANDIDATE / "genetics/figures/figure_manifest.tsv",
        pdf,
        source,
        CANDIDATE / "genetics/figures/sessionInfo.txt",
        out / "validation_results.tsv",
    ]
    seal = [{
        "release_id": RELEASE_ID, "gate": "GENETIC_READY", "status": "READY",
        "artifact": path.relative_to(CANDIDATE).as_posix(), "sha256": sha256(path),
    } for path in artifacts]
    atomic(CANDIDATE / "GENETIC_READY", ("release_id", "gate", "status", "artifact", "sha256"), seal)
    full = [{
        "release_id": RELEASE_ID,
        "gate": "FULL_READY",
        "status": "READY",
        "non_genetic_ready_sha256": sha256(CANDIDATE / "NON_GENETIC_READY"),
        "genetic_ready_sha256": sha256(CANDIDATE / "GENETIC_READY"),
        "candidate_audit_ready_sha256": sha256(CANDIDATE / "CANDIDATE_AUDIT_READY"),
    }]
    atomic(
        CANDIDATE / "FULL_READY",
        (
            "release_id", "gate", "status", "non_genetic_ready_sha256",
            "genetic_ready_sha256", "candidate_audit_ready_sha256",
        ),
        full,
    )
    print("GENETIC_READY\nFULL_READY")


if __name__ == "__main__":
    main()
