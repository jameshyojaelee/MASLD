#!/usr/bin/env python3
"""Freeze Plan 43 and correct the Tier-1/2 genetic locus provenance."""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from bridge_common import (
    CANDIDATE_ID,
    CANDIDATE_ROOT,
    CONFIG_ROOT,
    PROJECT_ROOT,
    atomic_write_text,
    read_tsv,
    sha256_file,
    stable_json_sha256,
    write_tsv,
)


EXPECTED_CLASSES = {
    "genetic_only": 413,
    "disease_state_only": 1227,
    "convergent": 34,
    "neither": 13257,
    "indeterminate_not_jointly_testable": 16298,
}
EXPECTED_TIER12_STUDIES = 35
EXPECTED_CORRECTED_LOCI = 326
EXPECTED_DRIVER_CORRECTIONS = 80
LOCUS_WINDOW_BP = 1_000_000


def finite_float(value: object, default: float = float("-inf")) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def parse_top_snp(value: str) -> tuple[str, int] | None:
    parts = str(value).strip().replace("_", ":").split(":")
    if len(parts) < 2:
        return None
    try:
        position = int(parts[1])
    except ValueError:
        return None
    return parts[0].removeprefix("chr"), position


def copy_and_verify_inputs() -> list[dict[str, object]]:
    manifest: list[dict[str, object]] = []
    for row in read_tsv(CONFIG_ROOT / "frozen_inputs.tsv"):
        source = PROJECT_ROOT / row["relative_path"]
        if not source.is_file():
            raise FileNotFoundError(source)
        observed = sha256_file(source)
        if observed != row["expected_sha256"]:
            raise RuntimeError(
                f"Frozen input drift for {row['input_id']}: {observed} != {row['expected_sha256']}"
            )
        destination = CANDIDATE_ROOT / "frozen_inputs" / f"{row['input_id']}__{source.name}"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        if sha256_file(destination) != observed:
            raise RuntimeError(f"Snapshot verification failed: {destination}")
        manifest.append(
            {
                "input_id": row["input_id"],
                "source_path": row["relative_path"],
                "snapshot_path": str(destination.relative_to(PROJECT_ROOT)),
                "sha256": observed,
                "size_bytes": source.stat().st_size,
                "role": row["role"],
            }
        )
    return manifest


def load_tier12_records(
    phenotype_rows: list[dict[str, str]],
) -> tuple[dict[str, dict[str, str]], dict[str, list[dict[str, str]]], list[dict[str, object]]]:
    study_meta = {
        row["study_name"]: row
        for row in phenotype_rows
        if row["tier"] in {"1", "2"} and row["placement"] == "main"
    }
    if len(study_meta) != EXPECTED_TIER12_STUDIES:
        raise RuntimeError(f"Expected 35 Tier-1/2 studies, observed {len(study_meta)}")
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    source_manifest: list[dict[str, object]] = []
    for study in sorted(study_meta):
        source = (
            PROJECT_ROOT
            / "GWAS/finemapping/results/susie_coloc"
            / study
            / "susie_coloc_combined.csv"
        )
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = CANDIDATE_ROOT / "frozen_inputs/tier12_coloc" / f"{study}.csv"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        observed = sha256_file(source)
        if sha256_file(destination) != observed:
            raise RuntimeError(f"Tier-1/2 snapshot verification failed: {study}")
        source_manifest.append(
            {
                "study_name": study,
                "tier": study_meta[study]["tier"],
                "trait": study_meta[study]["trait"],
                "direct_or_proxy": study_meta[study]["direct_or_proxy"],
                "source_path": str(source.relative_to(PROJECT_ROOT)),
                "snapshot_path": str(destination.relative_to(PROJECT_ROOT)),
                "sha256": observed,
                "size_bytes": source.stat().st_size,
            }
        )
        with source.open(newline="", encoding="utf-8") as handle:
            for record in csv.DictReader(handle):
                gene = record.get("gene", "").strip()
                if gene:
                    record = dict(record)
                    record["tier12_study"] = study
                    by_gene[gene].append(record)
    return study_meta, by_gene, source_manifest


def best_tier12_record(records: list[dict[str, str]]) -> dict[str, str]:
    with_susie = [r for r in records if finite_float(r.get("PP.H4.susie")) > float("-inf")]
    if not with_susie:
        return {}
    return sorted(
        with_susie,
        key=lambda row: (
            -finite_float(row.get("PP.H4.susie")),
            -finite_float(row.get("PP.H4.abf")),
            row["tier12_study"],
            row.get("top_snp", ""),
        ),
    )[0]


def build_corrected_registry(
    classes: list[dict[str, str]],
    study_meta: dict[str, dict[str, str]],
    by_gene: dict[str, list[dict[str, str]]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], int]:
    genetic = [row for row in classes if row["static_class"] == "genetic_only"]
    provisional: list[dict[str, object]] = []
    corrected_classes: list[dict[str, object]] = []
    mismatches = 0
    corrected_by_gene: dict[str, dict[str, object]] = {}
    for row in genetic:
        gene = row["gene_symbol"]
        best = best_tier12_record(by_gene.get(gene, []))
        if not best or finite_float(best.get("PP.H4.susie")) <= 0.5:
            raise RuntimeError(f"No qualifying Tier-1/2 SuSiE record for genetic-only gene {gene}")
        study = best["tier12_study"]
        meta = study_meta[study]
        parsed = parse_top_snp(best.get("top_snp", ""))
        if parsed is None:
            raise RuntimeError(f"Unresolved Tier-1/2 lead SNP for {gene} / {study}")
        tier1_support = sorted(
            {
                r["tier12_study"]
                for r in by_gene[gene]
                if study_meta[r["tier12_study"]]["tier"] == "1"
                and finite_float(r.get("PP.H4.susie")) > 0.5
            }
        )
        if study != row["driving_gwas"]:
            mismatches += 1
        record: dict[str, object] = {
            "gene_symbol": gene,
            "ensembl_id": row["ensembl_genetic"],
            "tier12_driving_gwas": study,
            "tier12_driving_trait": meta["trait"],
            "tier12_driving_tier": meta["tier"],
            "tier12_direct_or_proxy": meta["direct_or_proxy"],
            "tier12_top_snp": best.get("top_snp", ""),
            "chromosome": parsed[0],
            "position": parsed[1],
            "tier12_susie_pp4": best.get("PP.H4.susie", ""),
            "tier12_abf_pp4": best.get("PP.H4.abf", ""),
            "any_direct_tier1_support": str(bool(tier1_support)).lower(),
            "tier1_supporting_gwas": ";".join(tier1_support),
            "full_portfolio_driving_gwas": row["driving_gwas"],
            "full_portfolio_driving_trait": row["driving_trait"],
            "driver_corrected": str(study != row["driving_gwas"]).lower(),
            "locus_resolution": "tier12_lead_snp_1mb",
        }
        provisional.append(record)
        corrected_by_gene[gene] = record

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in provisional:
        grouped[(str(row["tier12_driving_gwas"]), str(row["chromosome"]))].append(row)
    for (study, chrom), rows in grouped.items():
        rows.sort(key=lambda item: (int(item["position"]), str(item["ensembl_id"])))
        component = 1
        previous: int | None = None
        for row in rows:
            position = int(row["position"])
            if previous is not None and position - previous > LOCUS_WINDOW_BP:
                component += 1
            row["tier12_locus_uid"] = f"{study}:chr{chrom}:component{component:04d}"
            previous = position

    by_locus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in provisional:
        by_locus[str(row["tier12_locus_uid"])].append(row)
    for rows in by_locus.values():
        selected = sorted(
            rows,
            key=lambda item: (
                -finite_float(item["tier12_susie_pp4"]),
                -finite_float(item["tier12_abf_pp4"]),
                str(item["ensembl_id"]),
            ),
        )[0]
        for item in rows:
            item["n_genetic_genes_in_locus"] = len(rows)
            item["is_representative"] = str(item is selected).lower()
            item["representative_gene"] = selected["gene_symbol"]

    for row in classes:
        output: dict[str, object] = dict(row)
        correction = corrected_by_gene.get(row["gene_symbol"])
        if correction:
            output.update(
                {
                    "tier12_driving_gwas": correction["tier12_driving_gwas"],
                    "tier12_driving_trait": correction["tier12_driving_trait"],
                    "tier12_driving_tier": correction["tier12_driving_tier"],
                    "tier12_top_snp": correction["tier12_top_snp"],
                    "tier12_locus_uid": correction["tier12_locus_uid"],
                    "any_direct_tier1_support": correction["any_direct_tier1_support"],
                    "full_portfolio_driving_gwas": correction["full_portfolio_driving_gwas"],
                    "full_portfolio_driving_trait": correction["full_portfolio_driving_trait"],
                    "driver_corrected": correction["driver_corrected"],
                }
            )
        corrected_classes.append(output)
    return sorted(provisional, key=lambda x: (str(x["tier12_locus_uid"]), str(x["ensembl_id"]))), corrected_classes, mismatches


def freeze_state_axis(registry: list[dict[str, str]]) -> list[dict[str, object]]:
    if len(registry) != 117:
        raise RuntimeError(f"Expected 117 programs, observed {len(registry)}")
    betas = [finite_float(row["primary_beta"], float("nan")) for row in registry]
    if any(not math.isfinite(beta) for beta in betas):
        raise RuntimeError("All 117 primary stage betas must be finite")
    denominator = sum(abs(beta) for beta in betas)
    if denominator <= 0:
        raise RuntimeError("Degenerate state-axis loading denominator")
    lineage_abs = Counter()
    for row, beta in zip(registry, betas):
        lineage_abs[row["cell_type"]] += abs(beta)
    output: list[dict[str, object]] = []
    for row, beta in zip(registry, betas):
        output.append(
            {
                "program_uid": row["program_uid"],
                "cell_type": row["cell_type"],
                "module": row["module"],
                "module_name": row["module_name"],
                "membership_sha256": row["membership_sha256"],
                "stage_beta": beta,
                "primary_loading": beta / denominator,
                "sign_loading": 1 if beta > 0 else -1 if beta < 0 else 0,
                "lineage_balanced_loading": (
                    beta / lineage_abs[row["cell_type"]] / len(lineage_abs)
                    if lineage_abs[row["cell_type"]] > 0
                    else 0
                ),
                "stage_model": row["model"],
                "stage_score_definition": row["score_definition"],
                "registry_release_id": row["release_id"],
            }
        )
    return output


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite existing candidate: {CANDIDATE_ROOT}")
    CANDIDATE_ROOT.mkdir(parents=True)
    input_manifest = copy_and_verify_inputs()
    classes = read_tsv(
        CANDIDATE_ROOT / "frozen_inputs/evidence_classes__frozen_evidence_classes.tsv"
    )
    counts = Counter(row["static_class"] for row in classes)
    if dict(counts) != EXPECTED_CLASSES:
        raise RuntimeError(f"Evidence-class drift: {dict(counts)}")
    phenotype = read_tsv(
        CANDIDATE_ROOT / "frozen_inputs/phenotype_registry__phenotype_registry.tsv"
    )
    study_meta, by_gene, coloc_manifest = load_tier12_records(phenotype)
    loci, corrected_classes, mismatches = build_corrected_registry(
        classes, study_meta, by_gene
    )
    n_loci = len({row["tier12_locus_uid"] for row in loci})
    if n_loci != EXPECTED_CORRECTED_LOCI:
        raise RuntimeError(f"Expected 326 corrected loci, observed {n_loci}")
    if mismatches != EXPECTED_DRIVER_CORRECTIONS:
        raise RuntimeError(f"Expected 80 driver corrections, observed {mismatches}")
    registry = read_tsv(
        CANDIDATE_ROOT / "frozen_inputs/program_registry__program_registry_v2.tsv"
    )
    axis = freeze_state_axis(registry)

    write_tsv(
        CANDIDATE_ROOT / "frozen_input_manifest.tsv",
        input_manifest,
        ["input_id", "source_path", "snapshot_path", "sha256", "size_bytes", "role"],
    )
    write_tsv(
        CANDIDATE_ROOT / "tier12_coloc_manifest.tsv",
        coloc_manifest,
        ["study_name", "tier", "trait", "direct_or_proxy", "source_path", "snapshot_path", "sha256", "size_bytes"],
    )
    locus_fields = list(loci[0])
    write_tsv(CANDIDATE_ROOT / "corrected_genetic_locus_registry.tsv", loci, locus_fields)
    class_fields = list(dict.fromkeys(field for row in corrected_classes for field in row))
    write_tsv(CANDIDATE_ROOT / "corrected_evidence_classes.tsv", corrected_classes, class_fields)
    write_tsv(
        CANDIDATE_ROOT / "frozen_state_axis.tsv",
        axis,
        list(axis[0]),
    )

    config_manifest: list[dict[str, object]] = []
    for source in sorted(CONFIG_ROOT.glob("*.tsv")):
        destination = CANDIDATE_ROOT / "frozen_spec" / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        config_manifest.append(
            {
                "source": str(source.relative_to(PROJECT_ROOT)),
                "snapshot": str(destination.relative_to(PROJECT_ROOT)),
                "sha256": sha256_file(source),
                "size_bytes": source.stat().st_size,
            }
        )
    write_tsv(
        CANDIDATE_ROOT / "config_manifest.tsv",
        config_manifest,
        ["source", "snapshot", "sha256", "size_bytes"],
    )

    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    specification = {
        "candidate_id": CANDIDATE_ID,
        "sealed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "git_head_at_seal": git_head,
        "ultimate_goal": "test whether the established MASLD transcriptome is a chronic, histologically reversible multicellular state and whether oriented regulatory-risk perturbations shift that frozen state geometry",
        "retrospective_disclosure": "the whole-registry reversal hypothesis arose after Plan 41 two-program outcomes; GSE106737 and Plan 41 contexts are discovery evidence",
        "outcome_unseen_external_tests": ["GSE83452", "GSE48452", "GSE281160_state_projection", "PXD052787_evidence_class_projection"],
        "class_counts": EXPECTED_CLASSES,
        "n_programs": 117,
        "n_tier12_studies": len(study_meta),
        "n_genetic_rows": len(loci),
        "n_corrected_loci": n_loci,
        "n_driver_corrections": mismatches,
        "locus_rule": "within Tier-1/2 driving GWAS and chromosome, lead SNPs separated by <=1,000,000 bp form connected components",
        "representative_rule": "highest Tier-1/2 SuSiE PP.H4, then Tier-1/2 ABF PP.H4, then Ensembl ID",
        "state_axis_rule": "all 117 frozen primary stage betas; L1-normalized signed loading; no significance selection",
        "minimum_programs_testable": 90,
        "minimum_loading_mass": 0.80,
        "minimum_cropseq_permutations": 100000,
        "source_outcomes_loaded_before_seal": False,
        "canonical_outputs_mutated": False,
        "frozen_input_manifest_sha256": sha256_file(CANDIDATE_ROOT / "frozen_input_manifest.tsv"),
        "tier12_coloc_manifest_sha256": sha256_file(CANDIDATE_ROOT / "tier12_coloc_manifest.tsv"),
        "corrected_locus_registry_sha256": sha256_file(CANDIDATE_ROOT / "corrected_genetic_locus_registry.tsv"),
        "corrected_evidence_classes_sha256": sha256_file(CANDIDATE_ROOT / "corrected_evidence_classes.tsv"),
        "frozen_state_axis_sha256": sha256_file(CANDIDATE_ROOT / "frozen_state_axis.tsv"),
        "config_manifest_sha256": sha256_file(CANDIDATE_ROOT / "config_manifest.tsv"),
    }
    specification["specification_sha256"] = stable_json_sha256(specification)
    atomic_write_text(
        CANDIDATE_ROOT / "SEALED.json",
        json.dumps(specification, indent=2, sort_keys=True) + "\n",
    )
    atomic_write_text(
        CANDIDATE_ROOT / "SPECIFICATION_SHA256",
        specification["specification_sha256"] + "\n",
    )
    print(json.dumps(specification, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
