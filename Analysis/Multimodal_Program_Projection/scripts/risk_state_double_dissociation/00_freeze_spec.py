#!/usr/bin/env python3
"""Freeze Plan 42 inputs, hypotheses, locus units, and predecessor identity."""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from risk_state_common import (
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
PRIMARY_MODULES = {("hepatocytes", "8"), ("hepatocytes", "20")}
LOCUS_WINDOW_BP = 1_000_000


def finite_float(value: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result if math.isfinite(result) else float("nan")


def coloc_records() -> dict[tuple[str, str], list[dict[str, str]]]:
    records: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    pattern = "GWAS/finemapping/results/susie_coloc/*/susie_coloc_combined.csv"
    for path in sorted(PROJECT_ROOT.glob(pattern)):
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                gene = row.get("gene", "").strip()
                gwas = row.get("gwas_name", "").strip()
                if gene and gwas:
                    records[(gene, gwas)].append(row)
    return records


def parse_top_snp(row: dict[str, str]) -> tuple[str, int] | None:
    top = row.get("top_snp", "").strip()
    if not top or top.upper() == "NA":
        return None
    parts = top.replace("_", ":").split(":")
    if len(parts) < 2:
        return None
    chrom = parts[0].removeprefix("chr")
    try:
        position = int(parts[1])
    except ValueError:
        return None
    return chrom, position


def build_locus_registry(classes: list[dict[str, str]]) -> list[dict[str, object]]:
    source = coloc_records()
    genetic = [row for row in classes if row["static_class"] == "genetic_only"]
    provisional: list[dict[str, object]] = []
    for row in genetic:
        gene = row["gene_symbol"]
        gwas = row["driving_gwas"]
        candidates = source.get((gene, gwas), [])
        candidates = sorted(
            candidates,
            key=lambda x: (
                -finite_float(x.get("PP.H4.susie", ""))
                if math.isfinite(finite_float(x.get("PP.H4.susie", "")))
                else float("inf"),
                -finite_float(x.get("PP.H4.abf", ""))
                if math.isfinite(finite_float(x.get("PP.H4.abf", "")))
                else float("inf"),
                x.get("top_snp", ""),
            ),
        )
        best = candidates[0] if candidates else {}
        parsed = parse_top_snp(best) if best else None
        provisional.append(
            {
                "gene_symbol": gene,
                "ensembl_id": row["ensembl_genetic"],
                "driving_gwas": gwas,
                "driving_trait": row["driving_trait"],
                "top_snp": best.get("top_snp", "") if best else "",
                "chromosome": parsed[0] if parsed else "",
                "position": parsed[1] if parsed else "",
                "source_susie_pp4": row["coloc_best_susie_pp4"],
                "source_abf_pp4": row["coloc_best_abf_pp4"],
                "locus_resolution": "lead_snp_1mb" if parsed else "unresolved_singleton",
            }
        )

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    unresolved: list[dict[str, object]] = []
    for row in provisional:
        if row["position"] == "":
            unresolved.append(row)
        else:
            grouped[(str(row["driving_gwas"]), str(row["chromosome"]))].append(row)

    for (gwas, chrom), rows in grouped.items():
        rows.sort(key=lambda x: (int(x["position"]), str(x["ensembl_id"])))
        component = 1
        previous: int | None = None
        for row in rows:
            position = int(row["position"])
            if previous is not None and position - previous > LOCUS_WINDOW_BP:
                component += 1
            row["locus_uid"] = f"{gwas}:chr{chrom}:component{component:04d}"
            previous = position
    for row in unresolved:
        row["locus_uid"] = f"{row['driving_gwas']}:unresolved:{row['ensembl_id']}"

    by_locus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in provisional:
        by_locus[str(row["locus_uid"])].append(row)
    for locus_rows in by_locus.values():
        def representative_key(item: dict[str, object]) -> tuple[float, float, str]:
            susie = finite_float(str(item["source_susie_pp4"]))
            abf = finite_float(str(item["source_abf_pp4"]))
            return (
                -(susie if math.isfinite(susie) else -1.0),
                -(abf if math.isfinite(abf) else -1.0),
                str(item["ensembl_id"]),
            )

        selected = sorted(locus_rows, key=representative_key)[0]
        for item in locus_rows:
            item["n_genetic_genes_in_locus"] = len(locus_rows)
            item["is_representative"] = str(item is selected).lower()
            item["representative_gene"] = selected["gene_symbol"]
    return sorted(provisional, key=lambda x: (str(x["locus_uid"]), str(x["ensembl_id"])))


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Candidate already exists; refusing overwrite: {CANDIDATE_ROOT}")
    CANDIDATE_ROOT.mkdir(parents=True)

    manifests: list[dict[str, object]] = []
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
            raise RuntimeError(f"Snapshot copy failed verification: {destination}")
        manifests.append(
            {
                "input_id": row["input_id"],
                "source_path": row["relative_path"],
                "snapshot_path": str(destination.relative_to(PROJECT_ROOT)),
                "sha256": observed,
                "size_bytes": source.stat().st_size,
                "role": row["role"],
            }
        )
    write_tsv(
        CANDIDATE_ROOT / "frozen_input_manifest.tsv",
        manifests,
        ["input_id", "source_path", "snapshot_path", "sha256", "size_bytes", "role"],
    )

    classes = read_tsv(CANDIDATE_ROOT / "frozen_inputs/evidence_classes__frozen_evidence_classes.tsv")
    counts = Counter(row["static_class"] for row in classes)
    if dict(counts) != EXPECTED_CLASSES:
        raise RuntimeError(f"Evidence-class counts drifted: {dict(counts)}")
    registry = read_tsv(CANDIDATE_ROOT / "frozen_inputs/program_registry__program_registry_v2.tsv")
    external = read_tsv(CANDIDATE_ROOT / "frozen_inputs/external_programs__external_test_programs.tsv")
    if len(registry) != 117:
        raise RuntimeError(f"Expected 117 frozen programs, observed {len(registry)}")
    if {(row["cell_type"], row["module"]) for row in external} != PRIMARY_MODULES:
        raise RuntimeError("External program family is not the frozen hepatocyte modules 8 and 20")

    loci = build_locus_registry(classes)
    if len(loci) != EXPECTED_CLASSES["genetic_only"]:
        raise RuntimeError(f"Expected 413 genetic-only locus rows, observed {len(loci)}")
    write_tsv(
        CANDIDATE_ROOT / "genetic_locus_registry.tsv",
        loci,
        [
            "locus_uid", "gene_symbol", "ensembl_id", "driving_gwas", "driving_trait",
            "top_snp", "chromosome", "position", "source_susie_pp4", "source_abf_pp4",
            "locus_resolution", "n_genetic_genes_in_locus", "is_representative",
            "representative_gene",
        ],
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

    amendment = {
        "amendment_type": "append_only_successor_source_identity",
        "predecessor": "public-functional-map-2026-08-09",
        "predecessor_mutated": False,
        "prior_skipped_object": {
            "expected_size_bytes": 4215858,
            "expected_md5": "fd3fcd19feec6572ff7df90c33f9d906",
            "status": "different_identity_skipped_by_plan41",
        },
        "official_final_object": {
            "source_id": "CLCC1_ST1",
            "expected_size_bytes": 4215920,
            "expected_md5": "c78f0e43640261bd77f26652b4fd2f93",
            "status": "source_gated_not_yet_loaded",
        },
        "prospective_status": "retrospective_execution_of_preexisting_plan41_question",
    }
    atomic_write_text(
        CANDIDATE_ROOT / "PLAN41_CLCC1_SOURCE_AMENDMENT.json",
        json.dumps(amendment, indent=2, sort_keys=True) + "\n",
    )

    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    specification = {
        "candidate_id": CANDIDATE_ID,
        "sealed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "git_head_at_seal": git_head,
        "class_counts": EXPECTED_CLASSES,
        "n_programs": len(registry),
        "n_external_programs": len(external),
        "n_genetic_rows": len(loci),
        "n_genetic_loci": len({row["locus_uid"] for row in loci}),
        "n_unresolved_singletons": sum(row["locus_resolution"] == "unresolved_singleton" for row in loci),
        "locus_rule": "within driving GWAS and chromosome, source lead colocalization SNPs separated by <=1,000,000 bp form connected components; unresolved genes are flagged singletons",
        "representative_rule": "highest source SuSiE PP.H4, then ABF PP.H4, then Ensembl ID",
        "minimum_interaction_permutations": 99999,
        "input_manifest_sha256": sha256_file(CANDIDATE_ROOT / "frozen_input_manifest.tsv"),
        "config_manifest_sha256": sha256_file(CANDIDATE_ROOT / "config_manifest.tsv"),
        "locus_registry_sha256": sha256_file(CANDIDATE_ROOT / "genetic_locus_registry.tsv"),
        "source_outcomes_loaded_before_seal": False,
        "canonical_outputs_mutated": False,
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
