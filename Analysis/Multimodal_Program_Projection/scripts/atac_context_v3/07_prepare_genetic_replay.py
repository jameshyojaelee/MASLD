#!/usr/bin/env python3
"""Fail-closed post-promotion COLOC replay planner and isolated runner builder."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

from atac_context_v3_lib import (
    RELEASE_ID,
    ContractError,
    candidate_root,
    default_candidate_root,
    project_root,
    sha256_file,
    write_tsv,
)


ROOT = project_root()
REQUIRED_MANIFEST = {
    "release_id",
    "status",
    "aggregate_path",
    "aggregate_sha256",
    "runner_path",
    "runner_sha256",
    "gwas_registry_path",
    "gwas_tier_path",
    "eqtl_susie_dir",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=default_candidate_root())
    parser.add_argument("--promoted-manifest", type=Path)
    parser.add_argument("--fixture-mode", action="store_true")
    return parser.parse_args()


def resolve(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def read_one_row(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not REQUIRED_MANIFEST.issubset(reader.fieldnames or []):
            raise ContractError(f"promoted manifest schema incomplete: {path}")
        rows = list(reader)
    if len(rows) != 1:
        raise ContractError("promoted COLOC manifest must contain exactly one row")
    return rows[0]


def require_hash(path: Path, expected: str) -> None:
    if not path.is_file() or sha256_file(path) != expected:
        raise ContractError(f"promoted input hash mismatch: {path}")


def write_blocked(genetics: Path, reason: str) -> None:
    genetics.mkdir(parents=True, exist_ok=True)
    write_tsv(
        genetics / "BLOCKED_UPSTREAM_RELEASE.tsv",
        ("release_id", "gate", "status", "reason", "created_utc"),
        [{
            "release_id": RELEASE_ID,
            "gate": "GENETIC_READY",
            "status": "blocked_upstream_release",
            "reason": reason,
            "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        }],
    )


def build_runner(source: Path, destination: Path) -> None:
    text = source.read_text(encoding="utf-8")
    out_old = 'OUT_DIR <- file.path(FM_DIR, paste0("results/susie_coloc", OUT_SUFFIX), gwas_name)'
    out_new = '''REPLAY_RUN_ROOT <- Sys.getenv("ATAC_V3_REPLAY_RUN_ROOT", unset = "")
if (REPLAY_RUN_ROOT == "") stop("ATAC_V3_REPLAY_RUN_ROOT is required")
OUT_DIR <- file.path(REPLAY_RUN_ROOT, gwas_name)'''
    filter_old = "egenes <- unique(eqtl_all$ENSG)"
    filter_new = '''egenes <- unique(eqtl_all$ENSG)
REPLAY_GENE_FILTER <- Sys.getenv("ATAC_V3_REPLAY_GENE_FILTER", unset = "")
if (REPLAY_GENE_FILTER == "" || !file.exists(REPLAY_GENE_FILTER)) {
  stop("ATAC_V3_REPLAY_GENE_FILTER is required")
}
replay_gene_ids <- unique(fread(REPLAY_GENE_FILTER)$ensembl)
egenes <- intersect(egenes, replay_gene_ids)
if (length(egenes) == 0L) stop("No requested replay genes occur on this chromosome")'''
    export_old = 'susie_method <- "susie"'
    export_new = '''susie_method <- "susie"
                replay_export <- Sys.getenv("ATAC_V3_REPLAY_EXPORT_ROOT", unset = "")
                if (replay_export == "") stop("ATAC_V3_REPLAY_EXPORT_ROOT is required")
                gene_export <- file.path(replay_export, gwas_name, paste0("chr", chr_num), gene_id)
                if (dir.exists(gene_export)) stop("Refusing to overwrite replay export: ", gene_export)
                dir.create(gene_export, recursive = TRUE, showWarnings = FALSE)
                signal_summary <- as.data.table(copy(susie_res$summary))
                signal_summary[, signal_pair_index := seq_len(.N)]
                signal_summary[, `:=`(
                  gwas_name = gwas_name, gene = gene_symbol, ensembl = gene_id,
                  chr = chr_num, gwas_signal = hit1, eqtl_signal = hit2
                )]
                posterior_wide <- as.data.table(copy(susie_res$results))
                posterior_cols <- setdiff(names(posterior_wide), "snp")
                posterior_long <- melt(
                  posterior_wide, id.vars = "snp", measure.vars = posterior_cols,
                  variable.name = "posterior_column", value.name = "SNP.PP.H4"
                )
                posterior_long[, signal_pair_index := match(posterior_column, posterior_cols)]
                allele_lookup <- unique(merged[, .(
                  snp = merge_key, hg19_position = eqtl_pos,
                  allele1 = gwas_a1, allele2 = gwas_a2
                )], by = "snp")
                posterior_long <- merge(posterior_long, allele_lookup, by = "snp", all.x = TRUE)
                posterior_long[, `:=`(
                  gwas_name = gwas_name, gene = gene_symbol, ensembl = gene_id,
                  chr = chr_num
                )]
                fwrite(signal_summary, file.path(gene_export, "signal_pairs.tsv"), sep = "\\t")
                fwrite(posterior_long, file.path(gene_export, "variant_posteriors.tsv.gz"), sep = "\\t")'''
    for old, new, label in (
        (out_old, out_new, "output isolation"),
        (filter_old, filter_new, "gene filter"),
        (export_old, export_new, "posterior export"),
    ):
        if text.count(old) != 1:
            raise ContractError(f"promoted runner cannot accept {label} patch exactly once")
        text = text.replace(old, new, 1)
    destination.parent.mkdir(parents=True, exist_ok=False)
    destination.write_text(text, encoding="utf-8")


def main() -> None:
    args = arguments()
    root = candidate_root(args.candidate_root, fixture_mode=args.fixture_mode)
    if not root.is_dir():
        raise ContractError("input freezer must create the candidate root first")
    genetics = root / "genetics"
    if args.promoted_manifest is None or not args.promoted_manifest.is_file():
        write_blocked(genetics, "promoted_corrected_coloc_manifest_unavailable")
        print("GENETIC_READY blocked: promoted corrected COLOC manifest unavailable")
        return

    manifest = read_one_row(args.promoted_manifest)
    if manifest["status"] != "PROMOTED":
        write_blocked(genetics, "corrected_coloc_manifest_not_promoted")
        print("GENETIC_READY blocked: corrected COLOC release is not promoted")
        return
    aggregate = resolve(manifest["aggregate_path"])
    runner = resolve(manifest["runner_path"])
    registry_path = resolve(manifest["gwas_registry_path"])
    tier_path = resolve(manifest["gwas_tier_path"])
    eqtl_dir = resolve(manifest["eqtl_susie_dir"])
    require_hash(aggregate, manifest["aggregate_sha256"])
    require_hash(runner, manifest["runner_sha256"])
    if not registry_path.is_file() or not tier_path.is_file() or not eqtl_dir.is_dir():
        raise ContractError("promoted COLOC manifest references missing frozen inputs")

    with tier_path.open("r", encoding="utf-8", newline="") as handle:
        tiers = list(csv.DictReader(handle, delimiter="\t"))
    primary = {
        row["study_name"]: (row["trait"], row["tier_label"])
        for row in tiers
        if row["placement"] == "main" and row["tier"] in {"1", "2"}
    }
    if len(primary) != 35:
        raise ContractError(f"expected 35 Tier-1/2 primary strata, observed {len(primary)}")
    with registry_path.open("r", encoding="utf-8", newline="") as handle:
        registry = {row["study_name"]: row for row in csv.DictReader(handle, delimiter="\t")}

    pairs = []
    with aggregate.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        needed = {"gwas_name", "gene", "ensembl", "chr", "PP.H4.susie", "method"}
        if not needed.issubset(reader.fieldnames or []):
            raise ContractError("promoted COLOC aggregate lacks replay columns")
        for row in reader:
            if row["gwas_name"] not in primary or row["method"] != "susie":
                continue
            try:
                pp4 = float(row["PP.H4.susie"])
            except (TypeError, ValueError):
                continue
            if pp4 <= 0.5:
                continue
            trait, trait_class = primary[row["gwas_name"]]
            pairs.append({
                "release_id": RELEASE_ID,
                "gwas_name": row["gwas_name"],
                "gene": row["gene"],
                "ensembl": row["ensembl"],
                "chr": int(row["chr"]),
                "trait": trait,
                "trait_class": trait_class,
                "ancestry": registry[row["gwas_name"]]["ancestry"],
                "promoted_pp_h4_susie": pp4,
            })
    pairs.sort(key=lambda row: (str(row["gwas_name"]), int(row["chr"]), str(row["ensembl"])))
    studies = sorted({str(row["gwas_name"]) for row in pairs})
    study_batch = {study: (index % 5) + 1 for index, study in enumerate(studies)}
    for row in pairs:
        row["batch_id"] = study_batch[str(row["gwas_name"])]
    columns = (
        "release_id", "batch_id", "gwas_name", "gene", "ensembl", "chr",
        "trait", "trait_class", "ancestry", "promoted_pp_h4_susie",
    )
    genetics.mkdir(parents=True, exist_ok=True)
    write_tsv(genetics / "replay_plan.tsv", columns, pairs)
    for batch in range(1, 6):
        write_tsv(
            genetics / "replay_batches" / f"batch_{batch}.tsv",
            columns,
            [row for row in pairs if row["batch_id"] == batch],
        )
    build_runner(runner, genetics / "replay_runner" / "06_susie_coloc_replay.R")
    write_tsv(
        genetics / "PROMOTED_UPSTREAM_GATE.tsv",
        ("release_id", "upstream_release_id", "status", "aggregate_sha256", "runner_sha256"),
        [{
            "release_id": RELEASE_ID,
            "upstream_release_id": manifest["release_id"],
            "status": "PROMOTED",
            "aggregate_sha256": manifest["aggregate_sha256"],
            "runner_sha256": manifest["runner_sha256"],
        }],
    )
    print(f"Prepared {len(pairs)} supported gene-study replay pairs in five deterministic batches")


if __name__ == "__main__":
    main()
