#!/usr/bin/env python3
"""Stage 4: MANIFEST.tsv for the runtime probe (inputs, code, outputs, with hashes)."""

from __future__ import annotations

import csv
import hashlib
import pathlib
import subprocess
import sys

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EXEC = PROJ / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z"
SRC = PROJ / "scripts/analysis/alphagenome_program"
WEIGHTS = PROJ / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-weights-20260915T103442Z"
FASTA = pathlib.Path(
    "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
    "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
)
HOSTED = (
    PROJ
    / "GWAS/finemapping/results/alphagenome_atlas/p6f-indel-rescue-20260914T135052Z/tables/indel_rescue_effects.tsv"
)


def sha256(path: pathlib.Path, limit_bytes: int | None = None) -> str:
    h = hashlib.sha256()
    n = 0
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
            n += len(chunk)
            if limit_bytes is not None and n >= limit_bytes:
                return h.hexdigest() + f"  (first {n} bytes only)"
    return h.hexdigest()


def rows():
    out = []

    def add(role, path, note=""):
        p = pathlib.Path(path)
        if not p.exists():
            out.append([role, str(p), "MISSING", "", note])
            return
        out.append([role, str(p), str(p.stat().st_size), sha256(p), note])

    add("weights_dir_inventory", WEIGHTS / "checkpoint_inventory.tsv",
        "12 files, 734,799,353 bytes, whole-artifact sha256 in acquisition.json")
    add("weights_acquisition", WEIGHTS / "acquisition.json",
        "google/alphagenome-all-folds @ a8f293a76ee73d5b57f3bf2ae146510589fcf187")
    for f in ("_CHECKPOINT_METADATA", "_METADATA", "manifest.ocdbt"):
        add("checkpoint_object", WEIGHTS / "checkpoints" / f, "Orbax OCDBT")
    add("reference_fasta", FASTA, "GRCh38 cellranger-arc-2024-A, hg38; the recipe's FASTA")
    add("reference_fasta_index", FASTA.with_suffix(".fa.fai"), "")
    add("hosted_reference_values", HOSTED,
        "hosted model-API run p6f-indel-rescue-20260914T135052Z; read-only, never modified")

    for p in sorted(SRC.glob("h1_*")):
        add("code", p, "written for this probe")

    code_repo = EXEC / "code/alphagenome_research"
    if code_repo.exists():
        rev = subprocess.run(
            ["git", "-C", str(code_repo), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=False).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(code_repo), "status", "--porcelain"],
            capture_output=True, text=True, check=False).stdout.strip()
        out.append(["code_repository", str(code_repo), "", rev,
                    f"google-deepmind/alphagenome_research pinned; worktree {'DIRTY' if dirty else 'clean'}"])

    for p in sorted(EXEC.glob("tables/*/*.json")):
        add("output", p, "")
    for p in sorted(EXEC.glob("env/*")):
        add("environment", p, "")
    return out


def main() -> None:
    dest = EXEC / "MANIFEST.tsv"
    r = rows()
    with dest.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["role", "path", "bytes", "sha256", "note"])
        w.writerows(r)
    print(f"wrote {dest} with {len(r)} rows")
    for row in r:
        print("\t".join(str(x) for x in row[:3]), file=sys.stderr)


if __name__ == "__main__":
    main()
