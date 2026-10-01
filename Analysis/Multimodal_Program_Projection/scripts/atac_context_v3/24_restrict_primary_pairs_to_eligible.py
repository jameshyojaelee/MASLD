#!/usr/bin/env python3
"""Restrict the Figure 5D primary-pair table to COLOC-eligible signal pairs.

Review item 2: the replay cut both SuSiE fits to the shared SNPs before
coloc.susie, so every signal pair looked eligible. The offline check
(GWAS/finemapping/src/perf/10_offline_eqtl_overlap_check.R) recomputed each
pair's shared posterior share on the full eQTL fit; the r2 COLOC release
reproduces those calls. This script applies them to the sealed replay context
without re-aggregating the 22 million variant posteriors.

09_aggregate_genetic_context.py chose the primary pair as the highest-PP.H4
pair over ALL pairs. Here it is chosen with the same key over ELIGIBLE pairs
only, so a gene-study whose best pair is now ineligible gets a different
primary pair, not just a filter. A gene-study leaves the table when no pair is
eligible (untestable) or its best eligible PP.H4 is <= 0.5. Per-pair lineage
masses are unchanged: for an eligible pair coloc restricts to the shared SNPs
itself, so the replay's SNP.PP.H4 vector is the one coloc uses.

Writes, under a NEW --out-root:
  genetics/context/genetic_lineage_context_primary_pairs.tsv  (09's schema)
  genetics/context/gene_study_eligibility.tsv
  genetics/context/restriction_summary.tsv
  genetics/context/restriction_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

from atac_context_v3_lib import (
    ContractError,
    LINEAGES,
    default_candidate_root,
    project_root,
    read_tsv,
    sha256_file,
    write_tsv,
)

OFFLINE = (
    project_root()
    / "GWAS/finemapping/results/coloc_eligibility_offline_20260923T211525Z"
)
CONTEXT = "genetics/context"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-root", type=Path, required=True,
                        help="new candidates/atac-context-v3-coloc-r2-<UTC ts>")
    parser.add_argument("--coloc-release", type=Path, required=True,
                        help="GWAS/finemapping/results/susie_coloc_r2_<ts> (ASSEMBLY_COMPLETE)")
    parser.add_argument("--eligibility", type=Path, default=OFFLINE / "eqtl_overlap_check.tsv")
    parser.add_argument("--gene-study", type=Path,
                        default=OFFLINE / "eqtl_overlap_check_gene_study.tsv")
    return parser.parse_args()


def primary_key(row: dict[str, str]) -> tuple:
    """09_aggregate_genetic_context.py: highest PP.H4, then signal ids, then index."""
    return (-float(row["pp_h4"]), row["gwas_signal"], row["eqtl_signal"],
            int(row["signal_pair_index"]))


def main() -> None:
    args = arguments()
    source_root = default_candidate_root()
    out_root = args.out_root.resolve()
    if out_root.parent != source_root.parent or not out_root.name.startswith(
        "atac-context-v3-coloc-r2-"
    ):
        raise ContractError(f"unsafe output root: {out_root}")
    if out_root.exists():
        raise ContractError(f"refusing to reuse {out_root}")
    if not (args.coloc_release / "ASSEMBLY_COMPLETE.json").is_file():
        raise ContractError("the COLOC release has not passed assembly")

    sealed = {row["artifact"]: row["sha256"]
              for row in read_tsv(source_root / "GENETIC_READY")[1]}
    inputs = {}
    for name in ("genetic_lineage_context_all_pairs.tsv",
                 "genetic_lineage_context_primary_pairs.tsv"):
        path = source_root / CONTEXT / name
        if sealed.get(f"{CONTEXT}/{name}") != sha256_file(path):
            raise ContractError(f"{name} does not match the sealed ATAC artifact")
        inputs[name] = path
    columns, all_pairs = read_tsv(inputs["genetic_lineage_context_all_pairs.tsv"])
    primary_columns, adopted_primary = read_tsv(inputs["genetic_lineage_context_primary_pairs.tsv"])
    if primary_columns != columns:
        raise ContractError("adopted primary and all-pair tables have different schemas")

    # (study, gene, idx1, idx2) -> signal_pair_index, from the replay exports.
    index_of: dict[tuple[str, str, str, str], int] = {}
    for path in sorted((source_root / "genetics/replay_execution").glob(
        "batch_*/exports/*/chr*/*/signal_pairs.tsv"
    )):
        for row in read_tsv(path)[1]:
            key = (row["gwas_name"], row["ensembl"], row["idx1"], row["idx2"])
            if key in index_of:
                raise ContractError(f"duplicated replay signal pair: {key}")
            index_of[key] = int(row["signal_pair_index"])

    by_pair: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    for row in all_pairs:
        by_pair[(row["gwas_name"], row["ensembl"], int(row["signal_pair_index"]))].append(row)
    if any(len(rows) != len(LINEAGES) for rows in by_pair.values()):
        raise ContractError("a signal pair does not have one row per lineage")

    eligible: dict[tuple[str, str, int], bool] = {}
    for row in read_tsv(args.eligibility)[1]:
        key = (row["gwas_name"], row["ensembl"], row["idx1"], row["idx2"])
        if key not in index_of:
            raise ContractError(f"eligibility row has no replay signal pair: {key}")
        pair = (key[0], key[1], index_of[key])
        if pair not in by_pair:
            raise ContractError(f"eligibility row has no lineage context: {pair}")
        if abs(float(row["replay_PP.H4"]) - float(by_pair[pair][0]["pp_h4"])) > 1e-12:
            raise ContractError(f"PP.H4 differs between the offline check and the context: {pair}")
        eligible[pair] = row["eligible"] == "TRUE"
    if set(eligible) != set(by_pair):
        raise ContractError("eligibility calls and lineage-context signal pairs differ")

    # The r2 release must make the same call on every pair it tested.
    with gzip.open(args.coloc_release / "susie_coloc_signal_pairs.tsv.gz", "rt", newline="") as handle:
        r2_calls = {
            (row["gwas"], row["ensembl"], row["idx1"], row["idx2"]): row["eligible"] == "TRUE"
            for row in csv.DictReader(handle, delimiter="\t")
        }
    disagree = [key for key, pair in ((k, (k[0], k[1], v)) for k, v in index_of.items())
                if pair in eligible and r2_calls.get(key) != eligible[pair]]
    if disagree:
        raise ContractError(f"{len(disagree)} pairs disagree with the r2 release, e.g. {disagree[:3]}")

    adopted_index = {(row["gwas_name"], row["ensembl"]): int(row["signal_pair_index"])
                     for row in adopted_primary}
    offline = {(row["gwas_name"], row["ensembl"]): row for row in read_tsv(args.gene_study)[1]}
    pairs_of: dict[tuple[str, str], list[int]] = defaultdict(list)
    for gwas, ensembl, index in by_pair:
        pairs_of[(gwas, ensembl)].append(index)
    if set(pairs_of) != set(offline) or set(pairs_of) != set(adopted_index):
        raise ContractError("gene-study pairs differ between context, offline check and primary table")

    primary_rows, gene_rows = [], []
    for key in sorted(pairs_of):
        lead = [by_pair[(key[0], key[1], index)][0] for index in pairs_of[key]
                if eligible[(key[0], key[1], index)]]
        best = min(lead, key=primary_key) if lead else None
        best_pp4 = float(best["pp_h4"]) if best else math.nan
        status = ("untestable_insufficient_shared_posterior" if best is None else
                  "retained_pp4_gt_0.5" if best_pp4 > 0.5 else "eligible_pp4_le_0.5")
        check = offline[key]
        if check["status"] != status or int(check["n_eligible"]) != len(lead) or \
                int(check["n_pairs"]) != len(pairs_of[key]) or \
                (best is not None and abs(float(check["corrected_pp4"]) - best_pp4) > 1e-12):
            raise ContractError(f"offline gene-study check disagrees for {key}")
        new_index = int(best["signal_pair_index"]) if best else ""
        old_rows = by_pair[(key[0], key[1], adopted_index[key])]
        gene_rows.append({
            "gwas_name": key[0], "gene": old_rows[0]["gene"], "ensembl": key[1],
            "trait_class": old_rows[0]["trait_class"], "n_pairs": len(pairs_of[key]),
            "n_eligible": len(lead), "adopted_primary_index": adopted_index[key],
            "adopted_pp_h4": old_rows[0]["pp_h4"], "r2_primary_index": new_index,
            "r2_pp_h4": best["pp_h4"] if best else "", "status": status,
            "primary_changed": str(best is not None and new_index != adopted_index[key]).upper(),
        })
        if status != "retained_pp4_gt_0.5":
            continue
        for row in by_pair[(key[0], key[1], new_index)]:
            primary_rows.append({**row, "primary_signal_pair": "TRUE"})

    # Unchanged primary pairs must reproduce the adopted rows exactly.
    adopted_rows = {(r["gwas_name"], r["ensembl"], r["lineage"]): r for r in adopted_primary}
    unchanged = {(r["gwas_name"], r["ensembl"]) for r in gene_rows
                 if r["status"] == "retained_pp4_gt_0.5" and r["primary_changed"] == "FALSE"}
    for row in primary_rows:
        if (row["gwas_name"], row["ensembl"]) in unchanged and \
                row != adopted_rows[(row["gwas_name"], row["ensembl"], row["lineage"])]:
            raise ContractError(f"unchanged primary pair does not reproduce: {row['gwas_name']} {row['ensembl']}")

    def census(rows: list[dict[str, str]]) -> dict[str, int]:
        pairs = {(r["gwas_name"], r["ensembl"], r["trait_class"]) for r in rows}
        out = {"lineage_rows": len(rows), "gene_study_pairs": len(pairs),
               "genes": len({r["ensembl"] for r in rows})}
        for trait_class in ("direct_MASLD", "liver_enzyme"):
            out[f"gene_study_pairs_{trait_class}"] = sum(p[2] == trait_class for p in pairs)
        return out

    old, new = census(adopted_primary), census(primary_rows)
    context = out_root / CONTEXT
    primary_path = context / "genetic_lineage_context_primary_pairs.tsv"
    write_tsv(primary_path, columns, primary_rows)
    write_tsv(context / "gene_study_eligibility.tsv", tuple(gene_rows[0]), gene_rows)
    status_counts = defaultdict(int)
    for row in gene_rows:
        status_counts[row["status"]] += 1
    summary = [{"metric": k, "adopted": old[k], "r2": new[k]} for k in old]
    summary += [{"metric": f"status_{k}", "adopted": "", "r2": v} for k, v in sorted(status_counts.items())]
    summary.append({"metric": "retained_primary_pair_changed", "adopted": "",
                    "r2": sum(r["primary_changed"] == "TRUE" and r["status"] == "retained_pp4_gt_0.5"
                              for r in gene_rows)})
    write_tsv(context / "restriction_summary.tsv", ("metric", "adopted", "r2"), summary)
    manifest = [{"role": role, "path": str(path), "sha256": sha256_file(path)} for role, path in (
        ("sealed_all_pairs", inputs["genetic_lineage_context_all_pairs.tsv"]),
        ("sealed_primary_pairs", inputs["genetic_lineage_context_primary_pairs.tsv"]),
        ("offline_signal_pairs", args.eligibility),
        ("offline_gene_study", args.gene_study),
        ("r2_signal_pairs", args.coloc_release / "susie_coloc_signal_pairs.tsv.gz"),
        ("restricted_primary_pairs", primary_path),
    )]
    write_tsv(context / "restriction_manifest.tsv", ("role", "path", "sha256"), manifest)
    print(json.dumps({"adopted": old, "r2": new, "status": dict(status_counts)}, indent=2))


if __name__ == "__main__":
    main()
