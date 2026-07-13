#!/usr/bin/env python3
"""Shared immutable contracts for SeqFunc v2 saturation mutagenesis."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
UNIFORM = Path(os.environ.get(
    "UNIFORM35_RUN_ROOT",
    ROOT / "GWAS/finemapping/runs/uniform35_v2_2026-07-13",
))
CHROMBPNET = Path(os.environ.get(
    "ADULT_CHROMBPNET_ROOT",
    ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2",
))
OUT = Path(os.environ.get(
    "SATURATION_V2_ROOT",
    ROOT / "GWAS/finemapping/results/seqfunc/haplotype_saturation/v2",
))

FOLD_GROUPS = {
    "fold0": ["chr1", "chr11", "chr15", "chr20", "chr21"],
    "fold1": ["chr2", "chr9", "chr14", "chr19", "chr22"],
    "fold2": ["chr3", "chr8", "chr10", "chr17"],
    "fold3": ["chr4", "chr7", "chr12", "chr18"],
    "fold4": ["chr5", "chr6", "chr13", "chr16"],
}
CHROM_TO_FOLD = {chrom: fold for fold, chroms in FOLD_GROUPS.items() for chrom in chroms}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def read_json(path: Path):
    with path.open() as handle:
        return json.load(handle)


def truthy(value) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "pass", "passed"}
