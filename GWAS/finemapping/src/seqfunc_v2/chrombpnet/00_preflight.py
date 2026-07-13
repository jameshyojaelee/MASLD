#!/usr/bin/env python3
"""Freeze donor, pseudoreplicate, reference, and five-fold contracts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from collections import defaultdict
from itertools import combinations
from pathlib import Path

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
SRC_RESULT = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet"
OUT = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"
RELEASE = ROOT / "GWAS/finemapping/results/seqfunc/releases/2026-07-13-r1"
ATAC = ROOT / "Analysis/ATAC/Human_Multiome"
ENV = ROOT / ".mamba/seqfunc_chrombpnet"
FASTA = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
BLACKLIST = Path("/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def executable(name: str) -> Path:
    path = ENV / "bin" / name
    if not path.is_file() or not os.access(path, os.X_OK):
        raise SystemExit(f"missing executable in isolated environment: {path}")
    return path


def read_counts(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["model_celltype"] == "hepatocyte"]
    if len(rows) != 18 or {r["donor_id"] for r in rows} != {f"D{i:02d}" for i in range(1, 19)}:
        raise SystemExit("expected exactly hepatocyte counts for D01-D18")
    return rows


def assign_pseudoreps(rows: list[dict[str, str]]) -> dict[str, int]:
    """Exactly balance donor counts and optimize cell mass within condition."""
    out: dict[str, int] = {}
    by_condition: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_condition[row["condition"]].append(row)
    if set(by_condition) != {"NORMAL", "MASL", "MASH"}:
        raise SystemExit(f"unexpected conditions: {sorted(by_condition)}")
    # The odd strata place their extra donor on opposite sides, producing exact
    # 9/9 global donor balance: NORMAL 3/2, MASL 2/2, MASH 4/5.
    rep1_targets = {"NORMAL": 3, "MASL": 2, "MASH": 4}
    for condition in sorted(by_condition):
        ordered = sorted(by_condition[condition], key=lambda x: x["donor_id"])
        total = sum(int(x["n_cells"]) for x in ordered)
        n_rep1 = rep1_targets[condition]
        candidates = []
        for idx in combinations(range(len(ordered)), n_rep1):
            chosen = {ordered[i]["donor_id"] for i in idx}
            rep1_cells = sum(int(ordered[i]["n_cells"]) for i in idx)
            candidates.append((abs(2 * rep1_cells - total), tuple(sorted(chosen)), chosen))
        _, _, best = min(candidates)
        for row in ordered:
            out[row["donor_id"]] = 1 if row["donor_id"] in best else 2
    if sum(x == 1 for x in out.values()) != 9 or sum(x == 2 for x in out.values()) != 9:
        raise SystemExit("pseudoreplicate assignment did not produce exact 9/9 donor balance")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(exist_ok=True)

    required = [
        RELEASE / "fold_manifest.tsv",
        RELEASE / "run_contract.json",
        SRC_RESULT / "cell_counts_by_donor.tsv",
        FASTA,
        Path(str(FASTA) + ".fai"),
        BLACKLIST,
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing/empty required input: {path}")
    for tool in ("chrombpnet", "macs3", "idr", "bedtools"):
        executable(tool)

    rows = read_counts(SRC_RESULT / "cell_counts_by_donor.tsv")
    assignment = assign_pseudoreps(rows)
    donor_rows = []
    for row in sorted(rows, key=lambda x: x["donor_id"]):
        donor = row["donor_id"]
        n_cells = int(row["n_cells"])
        barcode = SRC_RESULT / "barcodes" / f"{donor}.hepatocyte.txt"
        fragment = ATAC / "results/fragments" / f"{donor}_fragments.tsv.gz"
        if n_cells < 1000:
            raise SystemExit(f"donor {donor} has fewer than 1000 QC-passing hepatocytes")
        for path in (barcode, fragment):
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"missing donor input: {path}")
        n_barcodes = sum(1 for x in barcode.open() if x.strip())
        if n_barcodes != n_cells:
            raise SystemExit(f"{donor}: barcode count {n_barcodes} != cell count {n_cells}")
        donor_rows.append({
            "donor_id": donor,
            "condition": row["condition"],
            "n_hepatocytes": n_cells,
            "pseudoreplicate": assignment[donor],
            "barcode_file": str(barcode.resolve()),
            "fragment_file": str(fragment.resolve()),
        })

    manifest = out / "donor_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(donor_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(donor_rows)

    # Re-emit only autosomal sizes, which are the complete frozen model universe.
    sizes = {}
    with Path(str(FASTA) + ".fai").open() as handle:
        for line in handle:
            chrom, size = line.split("\t")[:2]
            if chrom in {f"chr{i}" for i in range(1, 23)}:
                sizes[chrom] = int(size)
    if set(sizes) != {f"chr{i}" for i in range(1, 23)}:
        raise SystemExit("reference index does not contain exactly chr1-chr22")
    with (out / "GRCh38.autosomes.chrom.sizes").open("w") as handle:
        for i in range(1, 23):
            handle.write(f"chr{i}\t{sizes[f'chr{i}']}\n")

    # Copy and validate the immutable fold table rather than recreating it.
    fold_src = RELEASE / "fold_manifest.tsv"
    fold_dst = out / "fold_manifest.tsv"
    shutil.copyfile(fold_src, fold_dst)
    with fold_dst.open() as handle:
        folds = list(csv.DictReader(handle, delimiter="\t"))
    expected_groups = {
        "0": "chr1,chr11,chr15,chr20,chr21",
        "1": "chr2,chr9,chr14,chr19,chr22",
        "2": "chr3,chr8,chr10,chr17",
        "3": "chr4,chr7,chr12,chr18",
        "4": "chr5,chr6,chr13,chr16",
    }
    if len(folds) != 5:
        raise SystemExit("release fold manifest does not contain five folds")
    for row in folds:
        if row["test_chromosomes"] != expected_groups[row["fold"]]:
            raise SystemExit(f"frozen fold {row['fold']} chromosome mismatch")

    versions = {}
    for tool, version_args in {
        "macs3": ["--version"],
        "idr": ["--version"],
        "bedtools": ["--version"],
    }.items():
        result = subprocess.run([str(executable(tool)), *version_args], text=True, capture_output=True, check=True)
        versions[tool] = (result.stdout or result.stderr).strip().splitlines()[0]
    versions["chrombpnet"] = "1.0.1 (environment lock; CLI has no version flag)"

    contract = {
        "status": "preflight_pass",
        "scope": "apply_only_supplementary",
        "model_cell_type": "adult_liver_hepatocyte",
        "n_donors": 18,
        "n_hepatocytes": sum(int(r["n_hepatocytes"]) for r in donor_rows),
        "conditions": {c: sum(r["condition"] == c for r in donor_rows) for c in ("NORMAL", "MASL", "MASH")},
        "pseudoreplicate_rule": "donor-intact exact 9/9 split; NORMAL 3/2, MASL 2/2, MASH 4/5; exhaustive within-condition minimum cell-mass imbalance; deterministic donor-id tie break",
        "peak_call": {"tool": "MACS3", "q": 0.05, "nomodel": True, "shift": -100, "extsize": 200, "input": "Tn5 insertion tagAlign"},
        "idr": {"inputs": "pseudoreplicates only", "rank": "signal.value", "global_idr_max": 0.05},
        "final_peak_rule": "pooled MACS3 peak retained iff its MACS3 summit falls inside a reproducible IDR interval",
        "chromosomes": [f"chr{i}" for i in range(1, 23)],
        "fold_source": str(fold_src.resolve()),
        "fold_sha256": sha256(fold_src),
        "reference_fasta": str(FASTA),
        "blacklist": str(BLACKLIST),
        "versions": versions,
        "seed": 42,
        "donor_heldout_claim_supported": False,
        "calibration_labels_used_for_training": False,
    }
    with (out / "run_contract.json").open("w") as handle:
        json.dump(contract, handle, indent=2, sort_keys=True)
        handle.write("\n")

    gate = {
        "fold0_internal_qc_pass": False,
        "currin_external_calibration_pass": False,
        "saturation_unlocked": False,
        "rule": "saturation_unlocked = fold0_internal_qc_pass AND currin_external_calibration_pass",
        "status": "awaiting_fold0_and_currin_evaluation",
    }
    gate_path = out / "gate_verdict.json"
    # Preflight is intentionally rerunnable, but must never relock or erase a
    # later reviewed verdict.
    if not gate_path.exists():
        with gate_path.open("w") as handle:
            json.dump(gate, handle, indent=2, sort_keys=True)
            handle.write("\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
