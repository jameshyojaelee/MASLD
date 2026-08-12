#!/usr/bin/env python3
"""Freeze ATAC v3 inputs and construct the donor-by-lineage QC contract."""

from __future__ import annotations

import argparse
import csv
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad

from atac_context_v3_lib import (
    COHORTS,
    DONOR_QC_COLUMNS,
    EXPECTED_PRIMARY_DONORS,
    INPUT_MANIFEST_COLUMNS,
    LABEL_MAP,
    LINEAGES,
    MIN_CELLS,
    PROGRAM_RELEASE_ID,
    RELEASE_ID,
    ContractError,
    candidate_root,
    default_candidate_root,
    project_root,
    sha256_file,
    sha256_text,
    write_json,
    write_tsv,
)


ROOT = project_root()
HOTSPOT = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / PROGRAM_RELEASE_ID
    / "hotspot"
)

COHORT_SPEC = {
    "GSE244832": {
        "h5ad": ROOT / "Analysis/ATAC/Human_Multiome/results/label_transfer/snapatac2_relabeled.h5ad",
        "metadata": ROOT / "data/GSE244832/metadata/donor_pairing.csv",
        "per_donor": ROOT / "Analysis/ATAC/Human_Multiome/results/snapatac2/per_donor",
        "fragments": ROOT / "Analysis/ATAC/Human_Multiome/results/fragments",
    },
    "GSE281367": {
        "h5ad": ROOT / "Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad",
        "metadata": ROOT / "data/GSE281367/metadata/donor_pairing.csv",
        "per_donor": ROOT / "Analysis/ATAC/Human_External/snapatac2_fast/per_donor",
        "fragments": ROOT / "Analysis/ATAC/Human_External/cellranger",
    },
}

REFERENCES = {
    "gencode_v49": Path("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"),
    "hg38_fasta": Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"),
    "hg38_fasta_index": Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa.fai"),
    "hg38_blacklist": Path("/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed"),
    "hg19_to_hg38_chain": ROOT / "data/broadaway_eqtl/hg19ToHg38.over.chain",
    "program_ready": HOTSPOT / "READY",
    "program_registry": HOTSPOT / "program_registry_v2.tsv",
    "program_membership": HOTSPOT / "program_membership_v2.tsv",
    "v2_dynamic_effects": ROOT / "Analysis/Multimodal_Program_Projection/results/atac/dynamic_program_accessibility.tsv",
    "v2_dynamic_scores": ROOT / "Analysis/Multimodal_Program_Projection/results/atac/dynamic_program_scores.tsv",
    "spatial_fig4f_source": ROOT / "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/figures/main/fig4_validation/data/fig4f_two_program_source_matrix.tsv",
}

EXISTING_INTEGRATION_TARGETS = (
    ROOT / "Analysis/ATAC/Integration/scripts/35_atac_integration.py",
    ROOT / "scripts/figures/fig4d_snatac_accessibility.R",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=default_candidate_root())
    parser.add_argument("--fixture-mode", action="store_true")
    parser.add_argument("--skip-large-hashes", action="store_true")
    return parser.parse_args()


def read_metadata(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle)
        return {str(row["donor_id"]): str(row["condition"]).upper() for row in rows}


def fragment_path(cohort: str, donor: str) -> Path:
    base = COHORT_SPEC[cohort]["fragments"]
    if cohort == "GSE244832":
        return base / f"{donor}_fragments.tsv.gz"
    return base / donor / "outs/fragments.tsv.gz"


def relative_or_absolute(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def file_row(path: Path, role: str, cohort: str = "", donor: str = "", skip_hash: bool = False) -> dict[str, object]:
    if not path.is_file():
        raise ContractError(f"missing input: {path}")
    return {
        "release_id": RELEASE_ID,
        "role": role,
        "cohort": cohort,
        "donor_id": donor,
        "relative_or_absolute_path": relative_or_absolute(path),
        "bytes": path.stat().st_size,
        "sha256": "deferred_large_file_hash" if skip_hash else sha256_file(path),
    }


def capture_preexisting(root: Path) -> None:
    rows = []
    diff_dir = root / "preexisting_changes"
    diff_dir.mkdir(parents=True, exist_ok=False)
    for index, path in enumerate(EXISTING_INTEGRATION_TARGETS, start=1):
        relative = path.relative_to(ROOT).as_posix()
        result = subprocess.run(
            ["git", "diff", "--", relative],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        diff_path = diff_dir / f"{index:02d}_{path.name}.diff"
        diff_path.write_text(result.stdout, encoding="utf-8")
        rows.append(
            {
                "path": relative,
                "sha256_before_v3": sha256_file(path),
                "preexisting_diff_sha256": sha256_text(result.stdout),
                "preexisting_diff_file": diff_path.relative_to(root).as_posix(),
            }
        )
    write_tsv(
        root / "preexisting_changes_manifest.tsv",
        ("path", "sha256_before_v3", "preexisting_diff_sha256", "preexisting_diff_file"),
        rows,
    )


def capture_environments(root: Path) -> None:
    destination = root / "environments"
    destination.mkdir(parents=True, exist_ok=False)
    micromamba = Path("/gpfs/commons/home/jameslee/.local/bin/micromamba")
    for environment in ("snapatac2", "rnaseq"):
        result = subprocess.run(
            [str(micromamba), "env", "export", "--explicit", "-n", environment],
            text=True,
            capture_output=True,
            check=True,
        )
        (destination / f"{environment}.explicit.txt").write_text(
            result.stdout, encoding="utf-8"
        )


def main() -> None:
    args = parse_args()
    out = candidate_root(args.candidate_root, fixture_mode=args.fixture_mode)
    if out.exists():
        raise ContractError(f"candidate already exists: {out}")
    out.mkdir(parents=True, exist_ok=False)

    manifest = []
    qc_rows = []
    package_root = Path(__file__).resolve().parent
    for path in sorted(package_root.rglob("*")):
        if (
            path.is_file()
            and "__pycache__" not in path.parts
            and path.suffix not in {".pyc", ".pyo"}
        ):
            manifest.append(file_row(path, "producer_source"))
    for role, path in REFERENCES.items():
        manifest.append(file_row(path, role))

    for cohort in COHORTS:
        spec = COHORT_SPEC[cohort]
        metadata = read_metadata(spec["metadata"])
        manifest.append(file_row(spec["h5ad"], "filtered_labeled_h5ad", cohort))
        manifest.append(file_row(spec["metadata"], "donor_metadata", cohort))
        data = ad.read_h5ad(spec["h5ad"], backed="r")
        obs = data.obs[["donor_id", "condition", "cell_type"]].copy()
        obs["lineage"] = obs["cell_type"].astype(str).map(LABEL_MAP[cohort])
        obs = obs.dropna(subset=["lineage"])
        grouped = obs.groupby(["donor_id", "lineage"], observed=True).size().to_dict()
        for donor, condition in metadata.items():
            per_donor = spec["per_donor"] / f"{donor}.h5ad"
            fragment = fragment_path(cohort, donor)
            manifest.append(file_row(per_donor, "fragment_backed_h5ad", cohort, donor))
            manifest.append(
                file_row(fragment, "raw_fragments", cohort, donor, args.skip_large_hashes)
            )
            manifest.append(
                file_row(Path(str(fragment) + ".tbi"), "raw_fragments_index", cohort, donor)
            )
            for lineage in LINEAGES:
                n_cells = int(grouped.get((donor, lineage), 0))
                eligible = n_cells >= MIN_CELLS
                qc_rows.append(
                    {
                        "release_id": RELEASE_ID,
                        "cohort": cohort,
                        "donor_id": donor,
                        "condition": condition,
                        "lineage": lineage,
                        "n_cells": n_cells,
                        "min_cells": MIN_CELLS,
                        "contrast_eligible": str(eligible).upper(),
                        "exclusion_reason": "" if eligible else "fewer_than_20_cells",
                    }
                )
        data.file.close()

    observed = {}
    for row in qc_rows:
        if row["contrast_eligible"] != "TRUE" or row["condition"] not in {"MASH", "NORMAL"}:
            continue
        key = (row["cohort"], row["lineage"])
        observed.setdefault(key, {"MASH": 0, "NORMAL": 0})[row["condition"]] += 1
    for key, expected in EXPECTED_PRIMARY_DONORS.items():
        actual = observed.get(key, {"MASH": 0, "NORMAL": 0})
        value = (actual["MASH"], actual["NORMAL"])
        if value != expected:
            raise ContractError(f"donor-lineage census drift for {key}: {value} != {expected}")

    write_tsv(out / "input_manifest.tsv", INPUT_MANIFEST_COLUMNS, manifest)
    write_tsv(out / "donor_lineage_qc.tsv", DONOR_QC_COLUMNS, qc_rows)
    write_tsv(
        out / "command_manifest.tsv",
        ("release_id", "stage", "command", "seed"),
        [
            {"release_id": RELEASE_ID, "stage": "freeze", "command": "python 01_freeze_inputs.py", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "peaks", "command": "python 02_call_consensus_peaks.py --n-jobs 16 --tempdir <scratch>", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "counts", "command": "python 03_count_fragments.py --cohort <cohort> --tempdir <scratch>", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "legacy_counts", "command": "python 03_count_fragments.py --cohort <cohort> --peak-space legacy_gse244832 --tempdir <scratch>", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "da", "command": "Rscript 04_run_da.R", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "programs", "command": "Rscript 05_run_programs.R", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "chromvar", "command": "Rscript 06_run_chromvar.R", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "condition_blindness", "command": "Rscript 10b_condition_blindness_audit.R", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "figures", "command": "Rscript 12_render_non_genetic.R", "seed": 42},
            {"release_id": RELEASE_ID, "stage": "validation", "command": "python 10a_validate_exact_recount.py && python 10_validate_non_genetic.py", "seed": 42},
        ],
    )
    capture_preexisting(out)
    capture_environments(out)
    write_json(
        out / "candidate_identity.json",
        {
            "release_id": RELEASE_ID,
            "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "seed": 42,
            "program_release_id": PROGRAM_RELEASE_ID,
            "external_outcomes_used_for_feature_definition": False,
            "canonical_release_modified": False,
        },
    )
    print(f"Wrote frozen ATAC v3 candidate inputs: {out}")


if __name__ == "__main__":
    main()
