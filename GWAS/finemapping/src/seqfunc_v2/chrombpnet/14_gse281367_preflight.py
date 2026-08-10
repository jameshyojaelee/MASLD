#!/usr/bin/env python3
"""Freeze the isolated GSE281367 hepatocyte ChromBPNet cohort contract."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import anndata as ad


ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
OUT = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2"
RELEASE = ROOT / "GWAS/finemapping/results/seqfunc/releases/2026-07-13-r1"
H5AD = ROOT / "Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad"
META = ROOT / "data/GSE281367/metadata/donor_pairing.csv"
PB_COL = ROOT / "Analysis/ATAC/Human_External/pseudobulk/hep_pseudobulk_coldata_GSE281367.tsv"
ENV = ROOT / ".mamba/seqfunc_chrombpnet"
FASTA = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
BLACKLIST = Path("/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed")
EXPECTED_DONORS = {f"Z{i:02d}" for i in range(1, 13)}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def executable(name: str) -> Path:
    path = ENV / "bin" / name
    if not path.is_file() or not os.access(path, os.X_OK):
        raise SystemExit(f"missing executable in isolated environment: {path}")
    return path


def original_fragment_barcode(obs_name: str) -> str:
    """Undo only AnnData's collision suffix, retaining the 10x '-1' suffix."""
    match = re.fullmatch(r"([ACGTN]+-1)(?:-\d+)?", obs_name)
    if match is None:
        raise SystemExit(f"unexpected GSE281367 cell identifier: {obs_name}")
    return match.group(1)


def assign_pseudoreps(rows: list[dict[str, object]]) -> dict[str, int]:
    """Split each 6-donor condition 3/3 while minimizing cell-mass imbalance."""
    assignment: dict[str, int] = {}
    by_condition: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_condition[str(row["condition"])].append(row)
    if set(by_condition) != {"NORMAL", "MASH"}:
        raise SystemExit(f"unexpected conditions: {sorted(by_condition)}")
    for condition in sorted(by_condition):
        ordered = sorted(by_condition[condition], key=lambda x: str(x["donor_id"]))
        if len(ordered) != 6:
            raise SystemExit(f"expected six {condition} donors, found {len(ordered)}")
        total = sum(int(x["n_hepatocytes"]) for x in ordered)
        candidates = []
        for idx in combinations(range(6), 3):
            chosen = {str(ordered[i]["donor_id"]) for i in idx}
            rep1_cells = sum(int(ordered[i]["n_hepatocytes"]) for i in idx)
            candidates.append((abs(2 * rep1_cells - total), tuple(sorted(chosen)), chosen))
        _, _, best = min(candidates)
        for row in ordered:
            donor = str(row["donor_id"])
            assignment[donor] = 1 if donor in best else 2
    if sum(x == 1 for x in assignment.values()) != 6 or sum(x == 2 for x in assignment.values()) != 6:
        raise SystemExit("pseudoreplicate assignment did not produce an exact 6/6 split")
    return assignment


def main() -> None:
    required = [H5AD, META, PB_COL, RELEASE / "fold_manifest.tsv", FASTA,
                Path(str(FASTA) + ".fai"), BLACKLIST]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing/empty required input: {path}")
    for tool in ("chrombpnet", "macs3", "idr", "bedtools"):
        executable(tool)

    with META.open() as handle:
        meta = list(csv.DictReader(handle))
    if {row["donor_id"] for row in meta} != EXPECTED_DONORS or len(meta) != 12:
        raise SystemExit("GSE281367 metadata does not contain exactly Z01-Z12")
    condition = {row["donor_id"]: row["condition"] for row in meta}
    if sum(x == "NORMAL" for x in condition.values()) != 6 or sum(x == "MASH" for x in condition.values()) != 6:
        raise SystemExit("GSE281367 condition contract is not 6 NORMAL / 6 MASH")

    atlas = ad.read_h5ad(H5AD, backed="r")
    for column in ("donor_id", "condition", "cell_type"):
        if column not in atlas.obs:
            raise SystemExit(f"source h5ad lacks obs.{column}")
    obs = atlas.obs.loc[atlas.obs["cell_type"].astype(str).eq("Hepatocyte"),
                        ["donor_id", "condition"]].copy()
    if len(obs) != 115_925 or set(obs["donor_id"].astype(str)) != EXPECTED_DONORS:
        raise SystemExit(f"expected 115925 hepatocyte labels across Z01-Z12, found {len(obs)}")
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["condition"] = obs["condition"].astype(str)
    for donor, group in obs.groupby("donor_id", sort=True):
        observed = set(group["condition"])
        if observed != {condition[donor]}:
            raise SystemExit(f"condition mismatch for {donor}: {observed} vs {condition[donor]}")

    with PB_COL.open() as handle:
        pb = {row["donor_id"]: int(row["n_cells"])
              for row in csv.DictReader(handle, delimiter="\t")}
    counts = obs.groupby("donor_id", sort=True).size().to_dict()
    if counts != pb:
        raise SystemExit("hepatocyte counts do not exactly reproduce the frozen pseudobulk coldata")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "logs").mkdir(exist_ok=True)
    barcode_dir = OUT / "barcodes"
    barcode_dir.mkdir(exist_ok=True)
    collision_suffixes_removed = 0
    donor_rows: list[dict[str, object]] = []
    for donor in sorted(EXPECTED_DONORS):
        names = obs.index[obs["donor_id"].eq(donor)].astype(str).tolist()
        barcodes = [original_fragment_barcode(x) for x in names]
        collision_suffixes_removed += sum(x != y for x, y in zip(names, barcodes))
        if len(set(barcodes)) != len(barcodes):
            raise SystemExit(f"normalizing AnnData collision suffixes created duplicate barcodes for {donor}")
        barcode_path = barcode_dir / f"{donor}.hepatocyte.txt"
        with barcode_path.open("w") as handle:
            handle.write("\n".join(sorted(barcodes)) + "\n")
        fragment = ROOT / f"Analysis/ATAC/Human_External/cellranger/{donor}/outs/fragments.tsv.gz"
        fragment_index = Path(str(fragment) + ".tbi")
        for path in (fragment, fragment_index):
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"missing donor fragment input: {path}")
        donor_rows.append({
            "donor_id": donor,
            "condition": condition[donor],
            "n_hepatocytes": counts[donor],
            "pseudoreplicate": 0,
            "barcode_file": str(barcode_path.resolve()),
            "fragment_file": str(fragment.resolve()),
            "fragment_bytes": fragment.stat().st_size,
        })
    assignment = assign_pseudoreps(donor_rows)
    for row in donor_rows:
        row["pseudoreplicate"] = assignment[str(row["donor_id"])]

    manifest = OUT / "donor_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(donor_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(donor_rows)

    sizes = {}
    with Path(str(FASTA) + ".fai").open() as handle:
        for line in handle:
            chrom, size = line.split("\t")[:2]
            if chrom in {f"chr{i}" for i in range(1, 23)}:
                sizes[chrom] = int(size)
    if set(sizes) != {f"chr{i}" for i in range(1, 23)}:
        raise SystemExit("reference index does not contain exactly chr1-chr22")
    with (OUT / "GRCh38.autosomes.chrom.sizes").open("w") as handle:
        for i in range(1, 23):
            handle.write(f"chr{i}\t{sizes[f'chr{i}']}\n")
    shutil.copyfile(RELEASE / "fold_manifest.tsv", OUT / "fold_manifest.tsv")

    versions = {}
    for tool, args in {"macs3": ["--version"], "idr": ["--version"],
                       "bedtools": ["--version"]}.items():
        result = subprocess.run([str(executable(tool)), *args], text=True,
                                capture_output=True, check=True)
        versions[tool] = (result.stdout or result.stderr).strip().splitlines()[0]
    versions["chrombpnet"] = "1.0.1 (environment lock; CLI has no version flag)"

    contract = {
        "status": "preflight_pass",
        "scope": "apply_only_cohort_comparison",
        "cohort": "GSE281367",
        "assay": "standalone_snATAC",
        "model_cell_type": "hepatocyte_gene_activity_marker_argmax",
        "cell_label_source": str(H5AD.resolve()),
        "cell_label_method": "04c positional repair; marker set and per-cell argmax copied from GSE244832 SnapATAC2 step5",
        "cell_label_rna_transfer_confirmed": False,
        "n_donors": 12,
        "n_hepatocytes": len(obs),
        "conditions": {name: sum(x == name for x in condition.values()) for name in ("NORMAL", "MASH")},
        "barcode_collision_suffixes_removed": collision_suffixes_removed,
        "pseudoreplicate_rule": "donor-intact 6/6 split; each replicate has 3 NORMAL and 3 MASH donors; exhaustive within-condition minimum cell-mass imbalance; deterministic donor-id tie break",
        "peak_call": {"tool": "MACS3", "q": 0.05, "nomodel": True, "shift": -100,
                      "extsize": 200, "input": "Tn5 insertion tagAlign"},
        "idr": {"inputs": "pseudoreplicates only", "rank": "signal.value", "global_idr_max": 0.05},
        "final_peak_rule": "pooled MACS3 peak retained iff its MACS3 summit falls inside a reproducible IDR interval",
        "chromosomes": [f"chr{i}" for i in range(1, 23)],
        "fold_source": str((RELEASE / "fold_manifest.tsv").resolve()),
        "fold_sha256": sha256(RELEASE / "fold_manifest.tsv"),
        "reference_fasta": str(FASTA),
        "blacklist": str(BLACKLIST),
        "source_metadata_sha256": sha256(META),
        "source_h5ad_bytes": H5AD.stat().st_size,
        "source_h5ad_mtime_ns": H5AD.stat().st_mtime_ns,
        "versions": versions,
        "seed": 42,
        "donor_heldout_claim_supported": False,
        "calibration_labels_used_for_training": False,
        "saturation_authority": False,
        "canonical_outputs_mutated": False,
    }
    with (OUT / "run_contract.json").open("w") as handle:
        json.dump(contract, handle, indent=2, sort_keys=True)
        handle.write("\n")
    gate_path = OUT / "gate_verdict.json"
    if not gate_path.exists():
        gate = {
            "all_internal_qc_pass": False,
            "currin_external_calibration_pass": False,
            "comparison_ready": False,
            "saturation_authority": False,
            "saturation_unlocked": False,
            "rule": "comparison_ready = all_internal_qc_pass AND currin_external_calibration_pass; saturation remains review-blocked",
            "status": "awaiting_training_and_currin_evaluation",
        }
        with gate_path.open("w") as handle:
            json.dump(gate, handle, indent=2, sort_keys=True)
            handle.write("\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
