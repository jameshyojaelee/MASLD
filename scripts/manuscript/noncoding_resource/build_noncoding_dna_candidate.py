#!/usr/bin/env python3
"""Build a pre-COLOC, credible-set-level noncoding-DNA candidate.

This producer deliberately excludes COLOC targets and evidence classes. It uses
the frozen July fine-mapping substrate to quantify variant consequence and
regulatory context while the corrected COLOC rerun is active. The final
synchronized genetics release must rebuild the same products from promoted
inputs before manuscript promotion.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from noncoding_genetics import (  # noqa: E402
    ABC_REQUIRED,
    ARCHITECTURE_FIELDS,
    ATAC_REQUIRED,
    CONSEQUENCE_REQUIRED,
    CS_REQUIRED,
    IDENTITY_REQUIRED,
    LIFTOVER_REQUIRED,
    MEMBER_FIELDS,
    REGULATORY_FIELDS,
    build_credible_set_products,
    build_identity_maps,
    gate_verdicts,
    load_trait_registry,
    load_unique_by_variant,
    parse_bool,
    read_table,
    require,
    sha256_file,
    write_tsv,
)

RELEASE_ID = "noncoding-dna-precoloc-candidate-2026-08-12"
SOURCE_ROLES = (
    "credible_sets",
    "consequence_annotation",
    "trait_registry",
    "gene_identity",
    "variant_liftover",
    "abc_context",
    "atac_context",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--credible-sets", type=Path, required=True)
    parser.add_argument("--consequence-annotation", type=Path, required=True)
    parser.add_argument("--trait-registry", type=Path, required=True)
    parser.add_argument("--gene-identity", type=Path, required=True)
    parser.add_argument("--variant-liftover", type=Path, required=True)
    parser.add_argument("--abc-context", type=Path, required=True)
    parser.add_argument("--atac-context", type=Path, required=True)
    return parser.parse_args()


def write_json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_source_manifest(path: Path, inputs: dict[str, Path]) -> None:
    rows = []
    for role in SOURCE_ROLES:
        source = inputs[role]
        require(source.is_file(), f"missing input: {source}")
        require(not source.is_symlink(), f"symlinked input is prohibited: {source}")
        rows.append(
            {
                "release_id": RELEASE_ID,
                "input_role": role,
                "source_path": str(source.resolve()),
                "size_bytes": str(source.stat().st_size),
                "sha256": sha256_file(source),
                "release_state": "frozen_candidate_not_promoted",
            }
        )
    write_tsv(
        path,
        rows,
        (
            "release_id",
            "input_role",
            "source_path",
            "size_bytes",
            "sha256",
            "release_state",
        ),
    )


def write_output_manifest(output_dir: Path) -> None:
    rows = []
    for path in sorted(output_dir.iterdir(), key=lambda item: item.name):
        if path.name == "output_manifest.tsv" or not path.is_file():
            continue
        rows.append(
            {
                "relative_path": path.name,
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256_file(path),
            }
        )
    write_tsv(
        output_dir / "output_manifest.tsv",
        rows,
        ("relative_path", "size_bytes", "sha256"),
    )


def main() -> None:
    args = parse_args()
    require(not args.output_dir.exists(), f"output exists: {args.output_dir}")
    inputs = {role: getattr(args, role) for role in SOURCE_ROLES}

    registry = load_trait_registry(args.trait_registry)
    identity_rows = read_table(args.gene_identity, IDENTITY_REQUIRED)
    _, _, tss_by_chromosome = build_identity_maps(identity_rows)
    consequences = load_unique_by_variant(
        args.consequence_annotation, CONSEQUENCE_REQUIRED
    )
    liftover = load_unique_by_variant(args.variant_liftover, LIFTOVER_REQUIRED)
    abc_enhancer_variants = {
        row["variant_id"]
        for row in read_table(args.abc_context, ABC_REQUIRED)
        if not parse_bool(
            row["abc_is_self_promoter"], f"ABC promoter {row['variant_id']}"
        )
    }
    atac_cell_types: dict[str, set[str]] = {}
    for row in read_table(args.atac_context, ATAC_REQUIRED):
        cell_type = row["cell_type"].strip()
        require(cell_type != "", f"empty ATAC cell type for {row['variant_id']}")
        atac_cell_types.setdefault(row["variant_id"], set()).add(cell_type)

    members, architecture, regulatory = build_credible_set_products(
        RELEASE_ID,
        read_table(args.credible_sets, CS_REQUIRED),
        registry,
        consequences,
        liftover,
        abc_enhancer_variants,
        atac_cell_types,
        tss_by_chromosome,
    )
    verdicts = gate_verdicts(RELEASE_ID, architecture)
    verdicts.extend(
        [
            {
                "release_id": RELEASE_ID,
                "record_type": "release_boundary",
                "trait_scope": "all_tier1_tier2",
                "method": "susie",
                "metric": "corrected_coloc_dependency",
                "value": "0",
                "threshold": "0",
                "passed": "true",
                "detail": "DNA architecture and accessibility do not use partial corrected COLOC outcomes",
            },
            {
                "release_id": RELEASE_ID,
                "record_type": "release_boundary",
                "trait_scope": "all_tier1_tier2",
                "method": "susie",
                "metric": "target_transcript_biotype_available",
                "value": "0",
                "threshold": "1",
                "passed": "false",
                "detail": "target-transcript biotype remains blocked until corrected COLOC promotion",
            },
        ]
    )

    temporary = args.output_dir.with_name(
        f".{args.output_dir.name}.tmp.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    )
    require(not temporary.exists(), f"temporary output exists: {temporary}")
    temporary.mkdir(parents=True, exist_ok=False)
    try:
        write_tsv(
            temporary / "credible_set_members_annotated.tsv", members, MEMBER_FIELDS
        )
        write_tsv(
            temporary / "credible_set_pip_architecture.tsv",
            architecture,
            ARCHITECTURE_FIELDS,
        )
        write_tsv(
            temporary / "credible_set_regulatory_context.tsv",
            regulatory,
            REGULATORY_FIELDS,
        )
        write_tsv(
            temporary / "genetics_noncoding_verdict.tsv",
            verdicts,
            (
                "release_id",
                "record_type",
                "trait_scope",
                "method",
                "metric",
                "value",
                "threshold",
                "passed",
                "detail",
            ),
        )
        write_source_manifest(temporary / "source_manifest.tsv", inputs)
        write_json(
            temporary / "execution_manifest.json",
            {
                "release_id": RELEASE_ID,
                "producer": str(Path(__file__).resolve()),
                "producer_sha256": sha256_file(Path(__file__).resolve()),
                "python_version": platform.python_version(),
                "n_reliable_credible_sets": len(architecture),
                "n_credible_set_members": len(members),
                "corrected_coloc_used": False,
                "target_gene_inference_included": False,
                "candidate_only": True,
            },
        )
        write_output_manifest(temporary)
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        temporary.rename(args.output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    print(
        json.dumps(
            {
                "output": str(args.output_dir),
                "n_credible_sets": len(architecture),
                "n_members": len(members),
                "corrected_coloc_used": False,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
