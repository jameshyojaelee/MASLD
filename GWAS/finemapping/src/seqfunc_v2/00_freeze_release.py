#!/usr/bin/env python3
"""Freeze the SeqFunc v2 release contract without touching canonical outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Iterable


RELEASE_ID = "2026-07-13-r1"
FOLD_GROUPS = {
    "G0": ["chr1", "chr11", "chr15", "chr20", "chr21"],
    "G1": ["chr2", "chr9", "chr14", "chr19", "chr22"],
    "G2": ["chr3", "chr8", "chr10", "chr17"],
    "G3": ["chr4", "chr7", "chr12", "chr18"],
    "G4": ["chr5", "chr6", "chr13", "chr16"],
}
FOLD_BP = {
    "G0": 597_188_383,
    "G1": 597_068_048,
    "G2": 560_489_058,
    "G3": 563_209_122,
    "G4": 557_046_911,
}


def sha256(path: Path, block: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(os.environ.get("MASLD_PROJECT_ROOT", Path(__file__).resolve().parents[4])),
    )
    parser.add_argument("--release-id", default=RELEASE_ID)
    args = parser.parse_args()
    root = args.root.resolve()
    fm = root / "GWAS/finemapping"
    out = fm / "results/seqfunc/releases" / args.release_id
    out.mkdir(parents=True, exist_ok=True)

    tier_path = fm / "config/gwas_trait_tier.tsv"
    registry_path = fm / "config/gwas_registry.tsv"
    tier = [r for r in read_tsv(tier_path) if r["placement"] == "main" and r["tier"] in {"1", "2"}]
    registry = {r["study_name"]: r for r in read_tsv(registry_path)}
    if len(tier) != 35 or {r["study_name"] for r in tier} - registry.keys():
        raise SystemExit("The frozen main Tier-1/2 portfolio is not exactly 35 registry-backed studies")

    studies = []
    for row in tier:
        reg = registry[row["study_name"]]
        ss = (fm / reg["sumstats_path"]).resolve()
        if not ss.is_file():
            raise SystemExit(f"Missing summary statistics: {ss}")
        ancestry = reg["ancestry"]
        block_scheme = {
            "EUR": "LDetect_EUR__PolyFun_UKBB_EUR_LD",
            "AFR": "LDetect_AFR__1KG_AFR_LD",
            "EAS": "LDetect_ASN__1KG_EAS_LD",
            "SAS": "LDetect_ASN__1KG_SAS_LD",
            "AMR": "LDetect_EUR_boundary_caveat__1KG_AMR_LD",
        }[ancestry]
        studies.append(
            {
                "study_name": row["study_name"],
                "trait": row["trait"],
                "tier": row["tier"],
                "tier_label": row["tier_label"],
                "ancestry": ancestry,
                "trait_type": reg["trait_type"],
                "N_tot": reg["N_tot"],
                "N_cases": reg["N_cases"],
                "sumstats_path": str(ss),
                "sumstats_size": ss.stat().st_size,
                "sumstats_mtime_ns": ss.stat().st_mtime_ns,
                "ld_panel": "polyfun_ukbb_eur" if ancestry == "EUR" else f"1kg_{ancestry.lower()}",
                "ld_panel_n": reg.get("ld_panel_n", ""),
                "block_scheme": block_scheme,
                "status": "frozen",
            }
        )
    studies.sort(key=lambda r: (int(r["tier"]), r["study_name"]))
    write_tsv(out / "study_manifest.tsv", studies, list(studies[0]))

    counts = {}
    for row in studies:
        counts[row["ancestry"]] = counts.get(row["ancestry"], 0) + 1
    if counts != {"EUR": 17, "EAS": 6, "AFR": 6, "SAS": 3, "AMR": 3}:
        raise SystemExit(f"Unexpected ancestry counts: {counts}")
    if sum(r["tier"] == "1" for r in studies) != 15:
        raise SystemExit("Unexpected Tier-1 count")

    folds = []
    names = list(FOLD_GROUPS)
    for i, test in enumerate(names):
        validation = names[(i + 1) % len(names)]
        train = [g for g in names if g not in {test, validation}]
        folds.append(
            {
                "fold": i,
                "seed": 42,
                "test_group": test,
                "test_chromosomes": ",".join(FOLD_GROUPS[test]),
                "test_bp": FOLD_BP[test],
                "validation_group": validation,
                "validation_chromosomes": ",".join(FOLD_GROUPS[validation]),
                "train_groups": ",".join(train),
                "train_chromosomes": ",".join(c for g in train for c in FOLD_GROUPS[g]),
            }
        )
    all_test = [chrom for row in folds for chrom in row["test_chromosomes"].split(",")]
    if sorted(all_test, key=lambda x: int(x.removeprefix("chr"))) != [f"chr{i}" for i in range(1, 23)]:
        raise SystemExit("Chromosome folds do not cover chr1-22 exactly once")
    write_tsv(out / "fold_manifest.tsv", folds, list(folds[0]))

    guarded = [
        root / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv",
        root / "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv",
        # The working checkout currently exposes the library build manifest and
        # composition summary but not the canonical CSV named in CLAUDE.md. Guard
        # the existing release authorities rather than fabricating the missing file.
        root / "Cas13_Library_Design/data/BUILD_MANIFEST.txt",
        root / "Cas13_Library_Design/results/library_composition_summary.csv",
        root / "RNA-seq/results/multi_evidence/convergence_evidence.csv",
    ]
    firewall = []
    for path in guarded:
        if not path.is_file():
            raise SystemExit(f"Canonical firewall target missing: {path}")
        firewall.append(
            {
                "path": str(path.resolve()),
                "size": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
                "sha256": sha256(path),
            }
        )
    write_tsv(out / "canonical_firewall_before.tsv", firewall, list(firewall[0]))

    claims = [
        {"claim_id": "uniform35_v2", "status": "predicted", "allowed_scope": "apply_only_sensitivity", "forbidden_scope": "canonical_coloc_or_atlas"},
        {"claim_id": "leafcutter_disease", "status": "predicted", "allowed_scope": "cluster_level_disease_association", "forbidden_scope": "inherited_splicing_without_sQTL"},
        {"claim_id": "hu_mpra", "status": "source_reproduced", "allowed_scope": "cell_line_reporter_assay", "forbidden_scope": "adult_hepatocyte_validation"},
        {"claim_id": "currin_caqtl", "status": "measured", "allowed_scope": "adult_bulk_liver_accessibility_QTL", "forbidden_scope": "primary_hepatocyte_or_expression_QTL"},
        {"claim_id": "chrombpnet_adult_liver", "status": "predicted", "allowed_scope": "chromosome_heldout_pooled_hepatocyte_model", "forbidden_scope": "donor_heldout_generalization"},
        {"claim_id": "saturation", "status": "nominated", "allowed_scope": "in_silico_mechanistic_nomination", "forbidden_scope": "endogenous_causal_effect"},
        {"claim_id": "abc_scenic_links", "status": "predicted", "allowed_scope": "source_separated_linked_PIP_mass", "forbidden_scope": "measured_link_or_gene_causal_probability"},
    ]
    write_tsv(out / "claim_ledger.tsv", claims, list(claims[0]))

    contract = {
        "release_id": args.release_id,
        "scope": "apply_only_supplementary",
        "uniform_prior": True,
        "functional_prior": False,
        "human_metafor": False,
        "canonical_promotion": False,
        "n_studies": 35,
        "tier_counts": {"1": 15, "2": 20},
        "ancestry_counts": counts,
        "chrombpnet_fold_type": "chromosome_heldout_not_donor_heldout",
        "second_seed_requires_user_review": True,
        "paths": {
            "finemap": str(fm / "runs/uniform35_v2_2026-07-13"),
            "splicing": str(fm / "results/seqfunc/disease_splicing/joint_v2"),
            "mpra": str(fm / "results/seqfunc/mpra_benchmark/v2"),
            "chrombpnet": str(fm / "results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"),
            "saturation": str(fm / "results/seqfunc/haplotype_saturation/v2"),
            "variant_to_gene": str(fm / "results/seqfunc/variant_to_gene/v2"),
        },
    }
    tmp = out / "run_contract.json.tmp"
    tmp.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, out / "run_contract.json")
    print(json.dumps({"release_root": str(out), "studies": len(studies), "folds": len(folds), "firewall_files": len(firewall)}, indent=2))


if __name__ == "__main__":
    main()
