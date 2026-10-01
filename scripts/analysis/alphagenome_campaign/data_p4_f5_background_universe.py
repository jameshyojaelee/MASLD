#!/usr/bin/env python3
"""Build the F5 null-background region universe, before any saturation score is read.

The original P4 universe (step 41) admitted an snATAC peak only when it carried a
disease-change state, cohort significance, genetic overlap or fixed117 program
linkage. F5's size-matched random gene-set null cannot draw from that universe:
it is selected, and a null drawn from a selected background is not a null. The
eligible background is every non-blacklisted consensus peak that the source DA
table names as a promoter peak of any gene, which is 24,626 regions over 17,834
genes; 18,570 of those were never registered and so were never retrieved.

This script re-derives that missing set from the two source tables and refuses to
continue unless it reproduces the independently checked deposit region for
region. Nothing here reads a saturation value, an effect or an outcome.

Output: <out>/tables/saturation_universe.tsv.gz in the schema step 42 consumes,
plus universe_receipt.json. Point AGA_OUT_ROOT at <out> to retrieve it.
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
ORIGINAL = ROOT / "GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z/tables/saturation_universe.tsv.gz"
DEPOSIT = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/p4-background-21775203/missing_background_regions.tsv.gz"
UNIVERSE_NAME = "U3_f5_background"
COLUMNS = ["universe", "region_key", "chrom", "start0", "end", "lineage", "u2a_da_state",
           "u2b_sig_either", "u2c_genetic_overlap", "u2d_program_linked", "width",
           "f5_background", "n_promoter_genes"]


def promoter_links():
    """gene -> set of peak coordinates, from the source DA table's promoter annotation only."""
    manifest = pd.read_csv(CTX / "consensus_peak_manifest.tsv", sep="\t", dtype=str)
    manifest = manifest.loc[~manifest.blacklist_overlap.str.upper().eq("TRUE")]
    eligible = set(manifest.lineage + "|" + manifest.chrom + ":" + manifest.start0 + "-" + manifest.end)
    source = pd.read_csv(CTX / "da/da_peak_results.tsv.gz", sep="\t", keep_default_na=False,
                         usecols=["lineage", "peak_coordinate", "promoter_genes"])
    source = source.loc[(source.lineage + "|" + source.peak_coordinate).isin(eligible)]
    links, lineages, genes_at = {}, {}, {}
    for row in source.itertuples(index=False):
        for gene in str(row.promoter_genes).replace(";", ",").split(","):
            if gene and gene not in {"NA", "nan", "none", "None"}:
                links.setdefault(gene, set()).add(row.peak_coordinate)
                lineages.setdefault(row.peak_coordinate, set()).add(row.lineage)
                genes_at.setdefault(row.peak_coordinate, set()).add(gene)
    return links, lineages, genes_at


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute-node execution required")
    args.out.mkdir(parents=True, exist_ok=False)
    links, lineages, genes_at = promoter_links()

    original = pd.read_csv(ORIGINAL, sep="\t", usecols=["universe", "region_key"])
    registered = set(original.loc[original.universe.eq("U2_snatac"), "region_key"])
    if len(registered) != 6930:
        raise ValueError("Original U2 coordinate universe changed")
    missing = set().union(*(r - registered for r in links.values()))

    # Exact-reproduction guard against the independently checked deposit.
    deposited = set(pd.read_csv(DEPOSIT, sep="\t").region_key)
    if missing != deposited:
        raise ValueError(f"Re-derivation differs from the checked deposit: "
                         f"{len(missing - deposited)} new, {len(deposited - missing)} absent")

    program_genes = set(pd.read_csv(PROGRAMS, sep="\t", usecols=["canonical_gene"]).canonical_gene.dropna())
    rows = []
    for key in sorted(missing):
        chrom, span = key.split(":")
        start0, end = (int(x) for x in span.split("-"))
        here = genes_at[key]
        rows.append({"universe": UNIVERSE_NAME, "region_key": key, "chrom": chrom, "start0": start0,
                     "end": end, "lineage": ";".join(sorted(lineages[key])), "u2a_da_state": "",
                     "u2b_sig_either": False, "u2c_genetic_overlap": False,
                     "u2d_program_linked": bool(here & program_genes), "width": end - start0,
                     "f5_background": True, "n_promoter_genes": len(here)})
    frame = pd.DataFrame(rows, columns=COLUMNS)
    if frame.region_key.duplicated().any() or len(frame) != len(missing):
        raise ValueError("Duplicate or dropped background coordinate")
    tables = args.out / "tables"
    tables.mkdir(parents=True)
    path = tables / "saturation_universe.tsv.gz"
    frame.to_csv(path, sep="\t", index=False)

    receipt = {
        "status": "f5_background_universe_defined_before_any_score_read",
        "universe": UNIVERSE_NAME,
        "regions": len(frame),
        "widths": sorted(frame.width.unique().tolist()),
        "total_bp": int(frame.width.sum()),
        "eligible_background_genes": len(links),
        "eligible_background_regions": len(set().union(*links.values())),
        "already_registered_excluded": len(registered),
        "program_linked_regions": int(frame.u2d_program_linked.sum()),
        "derivation": "source DA promoter_genes over non-blacklisted consensus peaks; no score, effect or outcome read",
        "reproduction_guard": "re-derived set equals the independently checked p4-background-21775203 deposit",
        "selection_note": "background is NOT restricted to registered or scored regions; that is the point of the family",
        "no_saturation_values_read": True,
        "universe_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (PROGRAMS, DEPOSIT)},
    }
    (args.out / "universe_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
