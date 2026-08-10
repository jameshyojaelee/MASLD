#!/usr/bin/env python3
"""Freeze separate Stage-A promoter and Stage-B exact-edit worklists.

Stage A models the genetically oriented target-expression change with a
promoter-directed gene perturbation: CRISPRa when risk increases expression and
CRISPRi when risk decreases it. Stage B edits an accessible shared-signal
variant exactly. This script performs reference-allele and sequence-window
checks only. It does not invent scores, install a design package, rank targets,
or freeze a target. Versioned CRISPick and PrimeDesign outputs remain required
external inputs for the actual guideability gate.
"""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    atomic_write_text,
    read_tsv,
    sha256_file,
    write_tsv,
)


FASTA = Path(
    "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
    "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
)
FAI = Path(f"{FASTA}.fai")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
WINDOW_RADIUS = 200


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def source_routing_root() -> Path:
    raw = os.environ.get("PLAN45_ROUTING_ROOT", "").strip()
    if not raw:
        raise RuntimeError("PLAN45_ROUTING_ROOT must point to a sealed RSR-03 routing release")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (
        PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
    ).resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"Routing source escapes candidate root: {path}")
    return path


class IndexedFasta:
    def __init__(self, fasta: Path, fai: Path):
        self.fasta = fasta
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with fai.open(encoding="utf-8") as handle:
            for raw in handle:
                name, length, offset, line_bases, line_bytes, *_ = raw.rstrip("\n").split("\t")
                self.index[name] = tuple(map(int, (length, offset, line_bases, line_bytes)))

    def fetch(self, chrom: str, start_1: int, end_1: int) -> str:
        if chrom not in self.index:
            raise RuntimeError(f"Chromosome absent from FASTA index: {chrom}")
        length, offset, line_bases, line_bytes = self.index[chrom]
        if start_1 < 1 or end_1 > length or end_1 < start_1:
            raise RuntimeError(f"Invalid FASTA interval: {chrom}:{start_1}-{end_1}")
        result = bytearray()
        cursor = start_1 - 1
        remaining = end_1 - start_1 + 1
        with self.fasta.open("rb") as handle:
            while remaining:
                line_offset = cursor % line_bases
                take = min(remaining, line_bases - line_offset)
                byte_offset = offset + (cursor // line_bases) * line_bytes + line_offset
                handle.seek(byte_offset)
                result.extend(handle.read(take))
                cursor += take
                remaining -= take
        sequence = result.decode("ascii").upper()
        if len(sequence) != end_1 - start_1 + 1 or set(sequence) - set("ACGTN"):
            raise RuntimeError(f"Invalid sequence retrieval: {chrom}:{start_1}-{end_1}")
        return sequence


def forward_alleles(
    reference: str, allele1: str, allele2: str, risk: str, palindromic: bool
) -> tuple[str, str, str]:
    if any(len(allele) != 1 or allele not in "ACGT" for allele in (allele1, allele2, risk)):
        return "", "", "failed_non_snv_or_invalid_allele"
    if risk not in {allele1, allele2} or allele1 == allele2:
        return "", "", "failed_risk_allele_not_biallelic"
    alleles = {allele1, allele2}
    if palindromic:
        return "", "", "unresolved_palindromic_forward_strand"
    if reference in alleles:
        forward_risk = risk
        forward_protective = allele2 if risk == allele1 else allele1
        return forward_risk, forward_protective, "direct_forward_match"
    complemented = {allele.translate(COMPLEMENT) for allele in alleles}
    if reference in complemented:
        forward_risk = risk.translate(COMPLEMENT)
        source_protective = allele2 if risk == allele1 else allele1
        return (
            forward_risk,
            source_protective.translate(COMPLEMENT),
            "complemented_to_forward_match",
        )
    return "", "", "failed_reference_allele_match"


def base_editor_route(reference: str, alternate: str) -> str:
    change = f"{reference}>{alternate}"
    if change in {"A>G", "T>C"}:
        return "ABE_candidate_reference_to_alternate"
    if change in {"C>T", "G>A"}:
        return "CBE_candidate_reference_to_alternate"
    return "prime_edit_required_reference_to_alternate"


def stage_a_mode(oriented_risk_effect: str) -> str:
    if oriented_risk_effect == "risk_increases_expression":
        return "CRISPRa"
    if oriented_risk_effect == "risk_decreases_expression":
        return "CRISPRi"
    raise RuntimeError(
        f"Cannot assign Stage-A perturbation mode from {oriented_risk_effect!r}"
    )


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite guide-design candidate: {CANDIDATE_ROOT}")
    routing_root = source_routing_root()
    seal_path = routing_root / "LINEAGE_ACCESSIBILITY_ROUTING_SEALED.json"
    route_path = routing_root / "lineage_accessibility_routing.tsv"
    variant_path = routing_root / "shared_signal_accessibility.tsv"
    for path in [seal_path, route_path, variant_path, FASTA, FAI]:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing guide-design source: {path}")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "lineage_accessibility_routing_complete_targets_not_selected":
        raise RuntimeError("Routing source is not a sealed target-freeze-prohibited release")
    if seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Routing source already claims frozen targets")
    if sha256_file(route_path) != seal["output_sha256"]["lineage_accessibility_routing"]:
        raise RuntimeError("Routing table hash mismatch")
    if sha256_file(variant_path) != seal["output_sha256"]["shared_signal_accessibility"]:
        raise RuntimeError("Routing variant hash mismatch")

    routes = {row["orientation_uid"]: row for row in read_tsv(route_path)}
    all_variants = read_tsv(variant_path)
    routed_uids = {uid for uid, row in routes.items() if yes(row["routing_gate_pass"])}
    exact_edit_candidates = [
        row
        for row in all_variants
        if row["orientation_uid"] in routed_uids and yes(row["accessibility_variant_pass"])
    ]
    if not exact_edit_candidates:
        raise RuntimeError("No routed variants are eligible for guide-design export")

    stage_a_rows: list[dict[str, object]] = []
    contracts: list[dict[str, object]] = []
    for index, uid in enumerate(sorted(routed_uids), start=1):
        route = routes[uid]
        mode = stage_a_mode(route["oriented_risk_effect"])
        request_uid = f"promoter-{index:04d}-{uid}"
        architecture_eligible = (
            route["dominant_celltype"] == "Hepatocytes"
            and route["accessibility_evidence_class"]
            == "two_cohort_hepatocyte_accessibility"
        )
        stage_a_rows.append(
            {
                "promoter_request_uid": request_uid,
                "orientation_uid": uid,
                "pair_family_uid": route["pair_family_uid"],
                "coarse_locus_uid": route["coarse_locus_uid"],
                "gene_symbol": route["gene_symbol"],
                "ensembl_id": route["ensembl_id"],
                "gwas_name": route["gwas_name"],
                "gwas_evidence_family": route["gwas_evidence_family"],
                "gwas_independence_class": route["gwas_independence_class"],
                "counts_for_replication_breadth": route[
                    "counts_for_replication_breadth"
                ],
                "trait": route["trait"],
                "tier": route["tier"],
                "phenotype_stratum": route["phenotype_stratum"],
                "susie_pp4": route["susie_pp4"],
                "orientation_consensus": route["orientation_consensus"],
                "dominant_celltype": route["dominant_celltype"],
                "accessibility_evidence_class": route[
                    "accessibility_evidence_class"
                ],
                "oriented_risk_effect": route["oriented_risk_effect"],
                "stage_a_perturbation_mode": mode,
                "crisprpick_gene_identifier": route["ensembl_id"]
                or route["gene_symbol"],
                "required_nonidentical_guides": 2,
                "primary_architecture_eligible_before_guide_gate": str(
                    architecture_eligible
                ).lower(),
                "guideability_status": "external_promoter_design_required",
                "target_freeze_status": (
                    "prohibited_until_external_design_and_protocol_gate"
                ),
            }
        )
        contracts.append(
            {
                "design_uid": request_uid,
                "orientation_uid": uid,
                "design_arm": "stage_a_promoter_gene_perturbation",
                "perturbation_mode": mode,
                "required_tool": "CRISPick",
                "required_tool_source": (
                    "https://portals.broadinstitute.org/gppx/crispick/public"
                ),
                "required_designs": 2,
                "input_field": "crisprpick_gene_identifier",
                "completion_rule": (
                    "at least two nonidentical source-tool-recommended promoter "
                    f"{mode} guides; full on/off-target export and tool version required"
                ),
            }
        )

    fasta = IndexedFasta(FASTA, FAI)
    stage_b_rows: list[dict[str, object]] = []
    fasta_records: list[str] = []
    for index, row in enumerate(
        sorted(
            exact_edit_candidates,
            key=lambda value: (
                value["orientation_uid"],
                -float(value["shared_posterior"]),
                value["snp"],
            ),
        ),
        start=1,
    ):
        chrom = row["chr_hg38"]
        position = int(row["pos_hg38"])
        start = position - WINDOW_RADIUS
        end = position + WINDOW_RADIUS
        sequence = fasta.fetch(chrom, start, end)
        reference = sequence[WINDOW_RADIUS]
        allele1 = row["gwas_effect_allele"].upper()
        allele2 = row["gwas_other_allele"].upper()
        risk = row["risk_allele"].upper()
        forward_risk, forward_protective, orientation_status = forward_alleles(
            reference,
            allele1,
            allele2,
            risk,
            yes(row["is_palindromic"]),
        )
        if orientation_status.startswith("failed") or orientation_status.startswith("unresolved"):
            alternate = ""
            installation_role = "unresolved"
            prime_input = ""
            editor_route = "unresolved"
            exact_edit_exportable = False
        else:
            if reference == forward_protective:
                alternate = forward_risk
                installation_role = "install_risk_from_protective_reference"
            elif reference == forward_risk:
                alternate = forward_protective
                installation_role = "install_protective_from_risk_reference"
            else:
                raise RuntimeError(f"Forward allele/reference inconsistency: {row['snp']}")
            prime_input = (
                sequence[:WINDOW_RADIUS]
                + f"({reference}/{alternate})"
                + sequence[WINDOW_RADIUS + 1 :]
            )
            editor_route = base_editor_route(reference, alternate)
            exact_edit_exportable = True

        design_uid = f"exact-{index:04d}-{row['orientation_uid']}-{chrom}-{position}"
        route = routes[row["orientation_uid"]]
        stage_b_rows.append(
            {
                "design_uid": design_uid,
                "orientation_uid": row["orientation_uid"],
                "pair_family_uid": route["pair_family_uid"],
                "coarse_locus_uid": route["coarse_locus_uid"],
                "gene_symbol": row["gene_symbol"],
                "ensembl_id": route["ensembl_id"],
                "gwas_name": route["gwas_name"],
                "gwas_evidence_family": route["gwas_evidence_family"],
                "gwas_independence_class": route["gwas_independence_class"],
                "counts_for_replication_breadth": route[
                    "counts_for_replication_breadth"
                ],
                "trait": route["trait"],
                "tier": route["tier"],
                "phenotype_stratum": route["phenotype_stratum"],
                "susie_pp4": route["susie_pp4"],
                "orientation_consensus": route["orientation_consensus"],
                "dominant_celltype": row["dominant_celltype"],
                "accessibility_evidence_class": route["accessibility_evidence_class"],
                "oriented_risk_effect": route["oriented_risk_effect"],
                "snp_hg19": row["snp"],
                "shared_posterior": row["shared_posterior"],
                "chromosome_hg38": chrom,
                "position_hg38": position,
                "window_start_hg38": start,
                "window_end_hg38": end,
                "reference_allele_forward": reference,
                "alternate_allele_forward": alternate,
                "risk_allele_forward": forward_risk,
                "protective_allele_forward": forward_protective,
                "forward_allele_orientation_status": orientation_status,
                "reference_to_alternate_role": installation_role,
                "exact_edit_exportable": str(exact_edit_exportable).lower(),
                "exact_edit_route": editor_route,
                "primedesign_input": prime_input,
                "sequence_401bp": sequence,
                "guideability_status": "external_design_required",
                "target_freeze_status": "prohibited_until_external_design_and_protocol_gate",
            }
        )
        fasta_records.append(f">{design_uid}|{chrom}:{start}-{end}\n{sequence}\n")
        contracts.append(
            {
                "design_uid": design_uid,
                "orientation_uid": row["orientation_uid"],
                "design_arm": "stage_b_exact_edit",
                "perturbation_mode": "exact_allele_edit",
                "required_tool": "PrimeDesign",
                "required_tool_source": "https://github.com/pinellolab/PrimeDesign",
                "required_designs": 2,
                "input_field": "primedesign_input",
                "completion_rule": (
                    "at least two nonidentical pegRNA designs for this exact variant with "
                    "complete PBS/RTT/ngRNA and off-target audit; a base-editor route is an "
                    "optional sensitivity"
                ),
            }
        )

    CANDIDATE_ROOT.mkdir(parents=True)
    stage_a_path = CANDIDATE_ROOT / "stage_a_promoter_design_worklist.tsv"
    stage_b_path = CANDIDATE_ROOT / "stage_b_exact_edit_worklist.tsv"
    contract_path = CANDIDATE_ROOT / "external_design_contract.tsv"
    fasta_path = CANDIDATE_ROOT / "stage_b_exact_edit_windows.fasta"
    write_tsv(stage_a_path, stage_a_rows, list(stage_a_rows[0]))
    write_tsv(stage_b_path, stage_b_rows, list(stage_b_rows[0]))
    write_tsv(contract_path, contracts, list(contracts[0]))
    atomic_write_text(fasta_path, "".join(fasta_records))
    manifest_rows = []
    for role, path in [
        ("routing_release", seal_path),
        ("routing_table", route_path),
        ("routing_variants", variant_path),
        ("GRCh38_reference_fasta", FASTA),
        ("GRCh38_reference_fasta_index", FAI),
    ]:
        manifest_rows.append(
            {
                "role": role,
                "source_path": (
                    str(path.relative_to(PROJECT_ROOT))
                    if PROJECT_ROOT in path.resolve().parents
                    else str(path)
                ),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    manifest_path = CANDIDATE_ROOT / "guide_design_input_manifest.tsv"
    write_tsv(manifest_path, manifest_rows, list(manifest_rows[0]))
    payload = {
        "status": "separate_stage_a_promoter_and_stage_b_exact_designs_pending",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_routing_root": str(routing_root.relative_to(PROJECT_ROOT)),
        "n_routed_pairs": len(routed_uids),
        "n_routed_pair_families": len(
            {routes[uid]["pair_family_uid"] for uid in routed_uids}
        ),
        "n_routed_physical_loci": len(
            {routes[uid]["coarse_locus_uid"] for uid in routed_uids}
        ),
        "n_stage_a_promoter_requests": len(stage_a_rows),
        "n_stage_b_accessible_variant_designs": len(stage_b_rows),
        "n_exact_edit_exportable": sum(
            yes(row["exact_edit_exportable"]) for row in stage_b_rows
        ),
        "sequence_window_bp": WINDOW_RADIUS * 2 + 1,
        "design_policy": (
            "one oriented promoter-level CRISPRa/i request per routed locus for Stage A; "
            "every routed accessibility-passing variant exported without ranking for "
            "Stage-B exact editing"
        ),
        "production_dependency_added": False,
        "scientific_outcomes_inspected": False,
        "experimental_targets_frozen": False,
        "next_gate": "external design results, off-target review, editability pilot, and protocol freeze",
        "output_sha256": {
            stage_a_path.name: sha256_file(stage_a_path),
            stage_b_path.name: sha256_file(stage_b_path),
            contract_path.name: sha256_file(contract_path),
            fasta_path.name: sha256_file(fasta_path),
            manifest_path.name: sha256_file(manifest_path),
        },
    }
    atomic_write_json(CANDIDATE_ROOT / "GUIDE_DESIGN_WORKLIST_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
