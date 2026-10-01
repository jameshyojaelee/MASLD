#!/usr/bin/env python3
"""Identify the original F5 missing null background without reading saturation.

Registered-region presence below is eligibility, not retrieval or integrity.
Gene and program membership are unchanged. This inventory never substitutes
the scored genes for the source's eligible random-gene-set universe.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CTX = ROOT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
PROGRAMS = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
UNIVERSE = ROOT / "GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z/tables/saturation_universe.tsv.gz"


def inventory(links, programs, registered):
    background = []
    for gene, regions in sorted(links.items()):
        background.append(dict(canonical_gene=gene, eligible_regions=len(regions),
            registered_regions=len(regions & registered), missing_regions=len(regions-registered)))
    rows = []
    for uid, part in programs.groupby("program_uid", sort=True):
        genes = set(part.canonical_gene)
        regions = set().union(*(links.get(g, set()) for g in genes))
        n = len(regions & registered)
        rows.append(dict(program_uid=uid, program_genes=len(genes), eligible_regions=len(regions),
            registered_regions=n, missing_regions=len(regions-registered),
            original_minimum200_possible_on_registered_universe=n >= 200,
            inference_status="metadata_only; no sensitivity or TF profile tested"))
    return pd.DataFrame(background), pd.DataFrame(rows)


def invariants():
    programs = pd.DataFrame(dict(program_uid=["fixed_program"], canonical_gene=["a"]))
    background, result = inventory({"a": {"chr1:10-20"}, "b": {"chr1:30-40"}}, programs, {"chr1:10-20"})
    assert len(background) == 2 and int(background.missing_regions.sum()) == 1
    assert result.iloc[0].registered_regions == 1
    assert not result.iloc[0].original_minimum200_possible_on_registered_universe
    return {"unscored_background_gene_preserved": True, "undersized_program_not_testable": True}


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute-node execution required")
    args.out.mkdir(parents=True, exist_ok=False)
    checks = invariants()
    manifest = pd.read_csv(CTX / "consensus_peak_manifest.tsv", sep="\t", dtype=str)
    manifest = manifest.loc[~manifest.blacklist_overlap.str.upper().eq("TRUE")]
    eligible = set(manifest.lineage + "|" + manifest.chrom + ":" + manifest.start0 + "-" + manifest.end)
    source = pd.read_csv(CTX / "da/da_peak_results.tsv.gz", sep="\t", keep_default_na=False,
                         usecols=["lineage", "peak_coordinate", "promoter_genes"])
    source = source.loc[(source.lineage + "|" + source.peak_coordinate).isin(eligible)]
    links = {}
    for row in source.itertuples(index=False):
        for gene in str(row.promoter_genes).replace(";", ",").split(","):
            if gene and gene not in {"NA", "nan", "none", "None"}:
                links.setdefault(gene, set()).add(row.peak_coordinate)
    programs = pd.read_csv(PROGRAMS, sep="\t", usecols=["program_uid", "canonical_gene"])
    if programs.program_uid.nunique() != 117 or programs.canonical_gene.isna().any():
        raise ValueError("Fixed117 source membership differs or has missing canonical genes")
    universe = pd.read_csv(UNIVERSE, sep="\t", usecols=["universe", "region_key"])
    registered = set(universe.loc[universe.universe.eq("U2_snatac"), "region_key"])
    if len(registered) != 6930:
        raise ValueError("Original U2 coordinate universe changed")
    background, result = inventory(links, programs, registered)
    missing = set().union(*(r-registered for r in links.values()))
    background.to_csv(args.out / "original_background_coverage.tsv", sep="\t", index=False)
    result.to_csv(args.out / "original_program_coverage.tsv", sep="\t", index=False)
    pd.DataFrame([dict(region_key=k, required_readout="same original hosted P4 saturation scorers")
                  for k in sorted(missing)]).to_csv(args.out / "missing_background_regions.tsv.gz", sep="\t", index=False)
    summary = dict(status="original_F5_information_requirement_quantified",
        eligible_background_genes=len(links), eligible_background_regions=len(set().union(*links.values())),
        missing_unique_background_regions=len(missing), background_genes_with_missing_regions=int(background.missing_regions.gt(0).sum()),
        programs=117, programs_reaching_original_minimum200=int(result.original_minimum200_possible_on_registered_universe.sum()),
        program_region_range=[int(result.registered_regions.min()), int(result.registered_regions.max())],
        no_saturation_values_read=True, retrieval_or_integrity_assessed=False, no_new_saturation_requested=True,
        background_restricted_to_registered_regions=False, scientific_F5_complete=False,
        fixed117_source_sha256=hashlib.sha256(PROGRAMS.read_bytes()).hexdigest(), scientific_invariants=checks,
        interpretation="Original random-gene-set null cannot use an outcome-selected or scored-only background; promoter annotation is restricted to the original source DA table")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
