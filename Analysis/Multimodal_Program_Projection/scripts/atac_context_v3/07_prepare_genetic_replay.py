#!/usr/bin/env python3
"""Fail-closed post-promotion COLOC replay planner and isolated runner builder."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
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
PACKAGE = Path(__file__).resolve().parent
REQUIRED_MANIFEST = {
    "release_id",
    "status",
    "aggregate_path",
    "aggregate_sha256",
    "runner_path",
    "runner_sha256",
    "gwas_registry_path",
    "gwas_registry_sha256",
    "gwas_tier_path",
    "gwas_tier_sha256",
    "eqtl_susie_dir",
    "upstream_input_manifest_path",
    "upstream_input_manifest_sha256",
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


def validate_upstream_input_manifest(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"path", "sha256"}.issubset(reader.fieldnames or []):
            raise ContractError("upstream input manifest requires path and sha256 columns")
        rows = list(reader)
    if not rows:
        raise ContractError("upstream input manifest is empty")
    for row in rows:
        require_hash(resolve(row["path"]), row["sha256"])
    return len(rows)


def write_blocked(genetics: Path, reason: str) -> None:
    genetics.mkdir(parents=True, exist_ok=True)
    path = genetics / "BLOCKED_UPSTREAM_RELEASE.tsv"
    if path.is_file():
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if (
            len(rows) == 1
            and rows[0].get("release_id") == RELEASE_ID
            and rows[0].get("gate") == "GENETIC_READY"
            and rows[0].get("status") == "blocked_upstream_release"
            and rows[0].get("reason") == reason
        ):
            return
        raise ContractError(f"refusing to replace an existing upstream gate: {path}")
    if path.is_symlink():
        raise ContractError(f"refusing to write through an upstream-gate symlink: {path}")
    write_tsv(
        path,
        ("release_id", "gate", "status", "reason", "created_utc"),
        [{
            "release_id": RELEASE_ID,
            "gate": "GENETIC_READY",
            "status": "blocked_upstream_release",
            "reason": reason,
            "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        }],
    )


def build_runner(source: Path, destination: Path, helper_source: Path | None = None) -> None:
    text = source.read_text(encoding="utf-8")
    helper_source = helper_source or PACKAGE / "genetics_export_helpers.R"
    if not helper_source.is_file():
        raise ContractError(f"missing replay export helper: {helper_source}")
    helper_old = "library(susieR)"
    helper_new = '''library(susieR)
REPLAY_HELPER <- Sys.getenv("ATAC_V3_REPLAY_HELPER", unset = "")
if (REPLAY_HELPER == "" || !file.exists(REPLAY_HELPER)) {
  stop("ATAC_V3_REPLAY_HELPER is required")
}
source(REPLAY_HELPER)'''
    registry_old = 'registry <- read.delim("config/gwas_registry.tsv", stringsAsFactors = FALSE)'
    registry_new = '''REPLAY_GWAS_REGISTRY <- Sys.getenv("ATAC_V3_GWAS_REGISTRY", unset = "")
if (REPLAY_GWAS_REGISTRY == "" || !file.exists(REPLAY_GWAS_REGISTRY)) {
  stop("ATAC_V3_GWAS_REGISTRY is required")
}
registry <- read.delim(REPLAY_GWAS_REGISTRY, stringsAsFactors = FALSE)'''
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
                export_coloc_susie_result(
                  susie_res = susie_res, merged = merged,
                  gwas_name = gwas_name, gene_symbol = gene_symbol,
                  gene_id = gene_id, chr_num = chr_num,
                  export_root = replay_export
                )'''
    for old, new, label in (
        (helper_old, helper_new, "export helper"),
        (registry_old, registry_new, "registry isolation"),
        (out_old, out_new, "output isolation"),
        (filter_old, filter_new, "gene filter"),
        (export_old, export_new, "posterior export"),
    ):
        if text.count(old) != 1:
            raise ContractError(f"promoted runner cannot accept {label} patch exactly once")
        text = text.replace(old, new, 1)
    destination.parent.mkdir(parents=True, exist_ok=False)
    destination.write_text(text, encoding="utf-8")
    shutil.copy2(helper_source, destination.parent / "genetics_export_helpers.R")


def main() -> None:
    args = arguments()
    root = candidate_root(args.candidate_root, fixture_mode=args.fixture_mode)
    if not root.is_dir():
        raise ContractError("input freezer must create the candidate root first")
    if not args.fixture_mode and not (root / "CANDIDATE_AUDIT_READY").is_file():
        raise ContractError("CANDIDATE_AUDIT_READY is required before genetics preparation")
    genetics = root / "genetics"
    if args.promoted_manifest is None or not args.promoted_manifest.is_file():
        write_blocked(genetics, "promoted_corrected_coloc_manifest_unavailable")
        print("GENETIC_READY blocked: promoted corrected COLOC manifest unavailable")
        return

    promoted_manifest = args.promoted_manifest.resolve()
    manifest = read_one_row(promoted_manifest)
    if manifest["status"] != "PROMOTED":
        write_blocked(genetics, "corrected_coloc_manifest_not_promoted")
        print("GENETIC_READY blocked: corrected COLOC release is not promoted")
        return
    aggregate = resolve(manifest["aggregate_path"])
    runner = resolve(manifest["runner_path"])
    registry_path = resolve(manifest["gwas_registry_path"])
    tier_path = resolve(manifest["gwas_tier_path"])
    eqtl_dir = resolve(manifest["eqtl_susie_dir"])
    upstream_input_manifest = resolve(manifest["upstream_input_manifest_path"])
    require_hash(aggregate, manifest["aggregate_sha256"])
    require_hash(runner, manifest["runner_sha256"])
    require_hash(registry_path, manifest["gwas_registry_sha256"])
    require_hash(tier_path, manifest["gwas_tier_sha256"])
    require_hash(upstream_input_manifest, manifest["upstream_input_manifest_sha256"])
    n_upstream_inputs = validate_upstream_input_manifest(upstream_input_manifest)
    if not eqtl_dir.is_dir():
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
    if {value[1] for value in primary.values()} != {"direct_MASLD", "liver_enzyme"}:
        raise ContractError("primary strata must contain both prespecified trait classes")
    with registry_path.open("r", encoding="utf-8", newline="") as handle:
        registry = {row["study_name"]: row for row in csv.DictReader(handle, delimiter="\t")}
    missing_registry = sorted(set(primary) - set(registry))
    if missing_registry or any(not registry[study].get("ancestry") for study in primary):
        raise ContractError(
            f"primary strata lack complete frozen registry metadata: {missing_registry[:5]}"
        )

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
    pair_keys = [(str(row["gwas_name"]), str(row["ensembl"])) for row in pairs]
    if not pairs or len(pair_keys) != len(set(pair_keys)):
        raise ContractError("primary SuSiE replay pairs are empty or duplicated")
    studies = sorted({str(row["gwas_name"]) for row in pairs})
    study_batch = {study: (index % 5) + 1 for index, study in enumerate(studies)}
    for row in pairs:
        row["batch_id"] = study_batch[str(row["gwas_name"])]
    columns = (
        "release_id", "batch_id", "gwas_name", "gene", "ensembl", "chr",
        "trait", "trait_class", "ancestry", "promoted_pp_h4_susie",
    )
    genetics.mkdir(parents=True, exist_ok=True)
    prepared = genetics / "prepared"
    if prepared.exists() or prepared.is_symlink():
        raise ContractError(f"refusing to overwrite prepared genetics stage: {prepared}")
    stale_pending = sorted(genetics.glob(".prepared.pending.*"))
    if stale_pending:
        raise ContractError(f"unresolved prior genetics preparation: {stale_pending}")
    pending = genetics / f".prepared.pending.{os.getpid()}"
    if pending.exists() or pending.is_symlink():
        raise ContractError(f"pending genetics preparation path exists: {pending}")
    pending.mkdir()
    write_tsv(pending / "replay_plan.tsv", columns, pairs)
    for batch in range(1, 6):
        write_tsv(
            pending / "replay_batches" / f"batch_{batch}.tsv",
            columns,
            [row for row in pairs if row["batch_id"] == batch],
        )
    build_runner(runner, pending / "replay_runner" / "06_susie_coloc_replay.R")
    direct_inputs = [
        ("promoted_terminal_manifest", promoted_manifest, promoted_manifest),
        ("promoted_aggregate", aggregate, aggregate),
        ("promoted_runner", runner, runner),
        ("gwas_registry", registry_path, registry_path),
        ("gwas_tier", tier_path, tier_path),
        ("upstream_input_manifest", upstream_input_manifest, upstream_input_manifest),
        (
            "isolated_replay_runner",
            pending / "replay_runner/06_susie_coloc_replay.R",
            prepared / "replay_runner/06_susie_coloc_replay.R",
        ),
        (
            "replay_export_helper",
            pending / "replay_runner/genetics_export_helpers.R",
            prepared / "replay_runner/genetics_export_helpers.R",
        ),
        ("replay_plan", pending / "replay_plan.tsv", prepared / "replay_plan.tsv"),
        *[
            (
                f"replay_batch_{batch}",
                pending / f"replay_batches/batch_{batch}.tsv",
                prepared / f"replay_batches/batch_{batch}.tsv",
            )
            for batch in range(1, 6)
        ],
        (
            "consensus_peak_manifest",
            root / "consensus_peak_manifest.tsv",
            root / "consensus_peak_manifest.tsv",
        ),
    ]
    requested_eqtl = sorted({
        eqtl_dir / f"chr{int(row['chr'])}" / f"{row['ensembl']}_susie.rds"
        for row in pairs
    })
    missing_eqtl = [str(path) for path in requested_eqtl if not path.is_file()]
    if missing_eqtl:
        raise ContractError(
            f"supported SuSiE pairs lack frozen eQTL objects: {missing_eqtl[:5]}"
        )
    input_rows = []
    for role, source_path, recorded_path in direct_inputs + [
        ("requested_eqtl_susie", path, path) for path in requested_eqtl
    ]:
        input_rows.append({
            "release_id": RELEASE_ID,
            "role": role,
            "path": str(recorded_path),
            "bytes": source_path.stat().st_size,
            "sha256": sha256_file(source_path),
        })
    write_tsv(
        pending / "replay_input_manifest.tsv",
        ("release_id", "role", "path", "bytes", "sha256"),
        input_rows,
    )
    prior_blocked = genetics / "BLOCKED_UPSTREAM_RELEASE.tsv"
    write_tsv(
        pending / "PROMOTED_UPSTREAM_GATE.tsv",
        (
            "release_id", "upstream_release_id", "status", "aggregate_path",
            "aggregate_sha256", "runner_path", "runner_sha256",
            "gwas_registry_path", "gwas_registry_sha256", "gwas_tier_path",
            "gwas_tier_sha256", "eqtl_susie_dir", "upstream_input_manifest_path",
            "upstream_input_manifest_sha256", "replay_plan_sha256",
            "replay_runner_sha256",
            "replay_input_manifest_sha256", "replay_helper_sha256",
            "n_upstream_inputs", "n_replay_pairs", "n_requested_eqtl_objects",
            "prior_blocked_gate_sha256",
        ),
        [{
            "release_id": RELEASE_ID,
            "upstream_release_id": manifest["release_id"],
            "status": "PROMOTED",
            "aggregate_path": str(aggregate),
            "aggregate_sha256": manifest["aggregate_sha256"],
            "runner_path": str(runner),
            "runner_sha256": manifest["runner_sha256"],
            "gwas_registry_path": str(registry_path),
            "gwas_registry_sha256": manifest["gwas_registry_sha256"],
            "gwas_tier_path": str(tier_path),
            "gwas_tier_sha256": manifest["gwas_tier_sha256"],
            "eqtl_susie_dir": str(eqtl_dir),
            "upstream_input_manifest_path": str(upstream_input_manifest),
            "upstream_input_manifest_sha256": manifest["upstream_input_manifest_sha256"],
            "replay_plan_sha256": sha256_file(pending / "replay_plan.tsv"),
            "replay_runner_sha256": sha256_file(
                pending / "replay_runner/06_susie_coloc_replay.R"
            ),
            "replay_input_manifest_sha256": sha256_file(pending / "replay_input_manifest.tsv"),
            "replay_helper_sha256": sha256_file(
                pending / "replay_runner/genetics_export_helpers.R"
            ),
            "n_upstream_inputs": n_upstream_inputs,
            "n_replay_pairs": len(pairs),
            "n_requested_eqtl_objects": len(requested_eqtl),
            "prior_blocked_gate_sha256": (
                sha256_file(prior_blocked) if prior_blocked.is_file() else ""
            ),
        }],
    )
    os.replace(pending, prepared)
    print(f"Prepared {len(pairs)} supported gene-study replay pairs in five deterministic batches")


if __name__ == "__main__":
    main()
