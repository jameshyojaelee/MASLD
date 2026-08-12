#!/usr/bin/env python3
"""Build noncoding genetic Resource tables from promoted inputs only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from noncoding_genetics import build_release


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--promotion-manifest", type=Path, required=True)
    parser.add_argument("--credible-sets", type=Path, required=True)
    parser.add_argument("--consequence-annotation", type=Path, required=True)
    parser.add_argument("--trait-registry", type=Path, required=True)
    parser.add_argument("--gene-identity", type=Path, required=True)
    parser.add_argument("--coloc", type=Path, required=True)
    parser.add_argument("--evidence-classes", type=Path, required=True)
    parser.add_argument("--variant-liftover", type=Path, required=True)
    parser.add_argument("--abc-context", type=Path, required=True)
    parser.add_argument("--atac-context", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = build_release(
        output_dir=args.output_dir,
        promotion_manifest=args.promotion_manifest,
        credible_sets=args.credible_sets,
        consequence_annotation=args.consequence_annotation,
        trait_registry=args.trait_registry,
        gene_identity=args.gene_identity,
        coloc=args.coloc,
        evidence_classes=args.evidence_classes,
        variant_liftover=args.variant_liftover,
        abc_context=args.abc_context,
        atac_context=args.atac_context,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
