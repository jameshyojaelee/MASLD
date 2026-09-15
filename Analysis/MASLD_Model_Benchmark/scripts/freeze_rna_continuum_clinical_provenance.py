"""Freeze and re-derive the clinical provenance of the five-cohort RNA continuum substrate.

The fibrosis and NAS values carried by ``five_cohort_sample_manifest.tsv`` for 290 of
its 844 samples (GSE130970 76, GSE135251 214) originate in two unchecksummed CSVs under
``RNA-seq/Human/Patient_Cohorts/archive/old_data_metadata/metadata/``.  Both cohorts are
``include_in_mega: true`` in ``config/human_datasets.yaml``, so losing that directory
hard-stops metadata harmonization at the BG-002 fail-closed check rather than degrading
quietly.  This script pins hashed copies inside the benchmark tree and re-derives every
clinical value from them so the join stops resting on an untraced assertion.

The unit of observation is the sequencing sample.  Four of the five cohorts carry no
donor key and none can be derived from the files on disk, so no output of this script
makes a donor-level claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

# The five cohorts of the continuum substrate, with the SraRunTable each one's clinical
# values are harmonized from.  Paths are relative to RNA-seq/Human/Patient_Cohorts and
# are transcribed from config/human_datasets.yaml `datasets.<GSE>.sra_table_path`.
COHORT_SOURCES: Mapping[str, str] = {
    "GSE126848": "pipelines/custom/GSE126848/metadata/SraRunTable.csv",
    "GSE130970": "archive/old_data_metadata/metadata/GSE130970_SraRunTable.csv",
    "GSE135251": "archive/old_data_metadata/metadata/GSE135251_SraRunTable.csv",
    "GSE162694": "pipelines/custom/GSE162694/metadata/SraRunTable.csv",
    "GSE213621": "pipelines/custom/GSE213621/metadata/SraRunTable.csv",
}

# Cohorts whose sra_table_path points into the archive directory.  These are the files
# the freeze exists to protect.
AT_RISK_COHORTS = ("GSE130970", "GSE135251")

# Fibrosis staging scale actually in force per cohort.  GSE213621 is NOT Kleiner 0-4:
# 00_harmonize_metadata.R:184-189 collapses Control->0, F0F1->1, F2->2, F3F4->3, so its
# level 1 pools F0 with F1, its level 3 pools F3 with F4, and it has no level 4 at all.
# Its top category is therefore a ceiling that aligns against the other cohorts'
# second-from-top level under a naive factor().  Never pool across the two scales.
FIBROSIS_SCALE_NATIVE: Mapping[str, str] = {
    "GSE126848": "not_recorded",
    "GSE130970": "kleiner_0_4",
    "GSE135251": "kleiner_0_4",
    "GSE162694": "kleiner_0_4_plus_normal_histology",
    "GSE213621": "gse213621_collapsed_0_3",
}

CROSS_COHORT_STAGE_COMPARABLE: Mapping[str, bool] = {
    "GSE126848": False,
    "GSE130970": True,
    "GSE135251": True,
    "GSE162694": True,
    "GSE213621": False,
}

# Provenance of the recorded-sex covariate.  Two cohorts covering 575 of 844 samples
# carry k-means-inferred sex rather than annotated sex, and both are the continuum's
# evaluation cohorts.  The frozen prespecification stratifies within-stage permutations
# on fibrosis_stage_by_inferred_sex and the primary model adjusts on it, so this field
# must travel with every downstream row.
SEX_SOURCE: Mapping[str, str] = {
    "GSE126848": "annotated",
    "GSE130970": "annotated",
    "GSE135251": "inferred_kmeans",
    "GSE162694": "annotated",
    "GSE213621": "inferred_kmeans",
}

NA_TOKENS = frozenset({"", "NA", "N/A", "na", "None", "null", "not applicable", "--"})


class ProvenanceError(RuntimeError):
    """A clinical value could not be re-derived from its frozen source."""


def is_missing(value: str | None) -> bool:
    return value is None or value.strip() in NA_TOKENS


def parse_optional_int(value: str | None) -> int | None:
    """Parse an integer clinical score, treating the recorded NA tokens as missing."""

    if is_missing(value):
        return None
    text = str(value).strip()
    try:
        # Some SRA exports carry integers as "3.0"; accept those without rounding
        # anything that is not already integral.
        number = float(text)
    except ValueError:
        return None
    if number != int(number):
        return None
    return int(number)


# --- Per-cohort re-derivation ------------------------------------------------------
#
# Each function mirrors one case of the harmonize() switch in
# analysis/integration/scripts/00_harmonize_metadata.R.  The line references name the
# block being reproduced so a reviewer can diff the two by eye.


def derive_gse126848(row: Mapping[str, str]) -> dict[str, Any]:
    """00_harmonize_metadata.R:97-113 -- no fibrosis or NAS column exists in source."""

    return {
        "fibrosis_stage": None,
        "nas_score": None,
        "fibrosis_native_field": "",
        "fibrosis_native_value": "",
        "nas_native_field": "",
        "nas_native_value": "",
    }


def derive_gse130970(row: Mapping[str, str]) -> dict[str, Any]:
    """00_harmonize_metadata.R:152-168."""

    return {
        "fibrosis_stage": parse_optional_int(row.get("fibrosis_stage")),
        "nas_score": parse_optional_int(row.get("nafld_activity_score")),
        "fibrosis_native_field": "fibrosis_stage",
        "fibrosis_native_value": (row.get("fibrosis_stage") or "").strip(),
        "nas_native_field": "nafld_activity_score",
        "nas_native_value": (row.get("nafld_activity_score") or "").strip(),
    }


def derive_gse135251(row: Mapping[str, str]) -> dict[str, Any]:
    """00_harmonize_metadata.R:131-149."""

    return {
        "fibrosis_stage": parse_optional_int(row.get("Fibrosis_stage")),
        "nas_score": parse_optional_int(row.get("nas_score")),
        "fibrosis_native_field": "Fibrosis_stage",
        "fibrosis_native_value": (row.get("Fibrosis_stage") or "").strip(),
        "nas_native_field": "nas_score",
        "nas_native_value": (row.get("nas_score") or "").strip(),
    }


def derive_gse162694(row: Mapping[str, str]) -> dict[str, Any]:
    """00_harmonize_metadata.R:196-219.

    "normal liver histology" maps to NA rather than to stage 0: normal histology is not
    the Kleiner F0 claim.  Those 31 samples are labelled Control and drop out of any
    factor(fibrosis_stage) fit, which is why this cohort contributes 109 and not 140.
    """

    native = (row.get("fibrosis_stage") or "").strip()
    if native == "normal liver histology":
        fibrosis: int | None = None
    else:
        fibrosis = parse_optional_int(native)
    return {
        "fibrosis_stage": fibrosis,
        "nas_score": parse_optional_int(row.get("nas_score")),
        "fibrosis_native_field": "fibrosis_stage",
        "fibrosis_native_value": native,
        "nas_native_field": "nas_score",
        "nas_native_value": (row.get("nas_score") or "").strip(),
    }


# 00_harmonize_metadata.R:184-189.  Note the absent level 4.
GSE213621_STAGE_RECODE: Mapping[str, int] = {
    "Control": 0,
    "F0F1": 1,
    "F2": 2,
    "F3F4": 3,
}


def derive_gse213621(row: Mapping[str, str]) -> dict[str, Any]:
    """00_harmonize_metadata.R:172-193 -- collapsed scale, and nas_score is NA cohort-wide."""

    native = (row.get("fibrotic_stage") or "").strip()
    return {
        "fibrosis_stage": GSE213621_STAGE_RECODE.get(native),
        "nas_score": None,
        "fibrosis_native_field": "fibrotic_stage",
        "fibrosis_native_value": native,
        "nas_native_field": "",
        "nas_native_value": "",
    }


DERIVERS: Mapping[str, Callable[[Mapping[str, str]], dict[str, Any]]] = {
    "GSE126848": derive_gse126848,
    "GSE130970": derive_gse130970,
    "GSE135251": derive_gse135251,
    "GSE162694": derive_gse162694,
    "GSE213621": derive_gse213621,
}


# --- Freezing ----------------------------------------------------------------------


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class FrozenInput:
    role: str
    source_path: Path
    frozen_path: Path
    sha256: str


def freeze_input(source: Path, destination: Path, role: str) -> FrozenInput:
    """Copy one read-only source file into the benchmark tree and record its digest."""

    if not source.is_file():
        raise ProvenanceError(f"clinical provenance source is missing: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    digest = sha256_file(destination)
    if digest != sha256_file(source):
        raise ProvenanceError(f"frozen copy does not match its source: {source}")
    return FrozenInput(role=role, source_path=source, frozen_path=destination, sha256=digest)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_cell(row.get(key)) for key in fieldnames})


def format_cell(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return str(value)


# --- Reproduction ------------------------------------------------------------------

REPRODUCTION_FIELDS = (
    "analysis_unit_id",
    "sample_id",
    "sra_run_accession",
    "dataset",
    "unit_of_observation",
    "n_donors",
    "n_donors_status",
    "fibrosis_manifest_value",
    "fibrosis_rederived_value",
    "fibrosis_agrees",
    "fibrosis_native_field",
    "fibrosis_native_value",
    "fibrosis_scale_native",
    "cross_cohort_stage_comparable",
    "nas_manifest_value",
    "nas_rederived_value",
    "nas_agrees",
    "nas_native_field",
    "nas_native_value",
    "sex_source",
    "source_table_role",
)


def build_reproduction(
    manifest_rows: Sequence[Mapping[str, str]],
    source_tables: Mapping[str, Sequence[Mapping[str, str]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Re-derive every manifest clinical value from the frozen source tables.

    Returns the per-sample reproduction rows and the rows present in a source table but
    absent from the manifest.  The second list is reported rather than explained: the
    yaml records ``excluded_samples: []`` for both archived cohorts, so these drops
    happen somewhere downstream of metadata harmonization and this script does not
    claim to know where.
    """

    indexed: dict[str, dict[str, Mapping[str, str]]] = {}
    for dataset, rows in source_tables.items():
        by_run: dict[str, Mapping[str, str]] = {}
        for row in rows:
            run = (row.get("Run") or "").strip()
            if not run:
                continue
            if run in by_run:
                raise ProvenanceError(
                    f"{dataset}: join key Run is not unique, found {run} twice"
                )
            by_run[run] = row
        indexed[dataset] = by_run

    reproduction: list[dict[str, Any]] = []
    for entry in manifest_rows:
        dataset = entry["dataset"]
        run = entry["sample_id"]
        source_row = indexed.get(dataset, {}).get(run)
        if source_row is None:
            raise ProvenanceError(
                f"{dataset}: manifest sample {run} has no row in the frozen source table"
            )
        derived = DERIVERS[dataset](source_row)
        manifest_fibrosis = parse_optional_int(entry.get("fibrosis_stage"))
        manifest_nas = parse_optional_int(entry.get("nas_score"))
        reproduction.append(
            {
                "analysis_unit_id": entry.get("analysis_unit_id", f"{dataset}::{run}"),
                "sample_id": run,
                "sra_run_accession": run,
                "dataset": dataset,
                "unit_of_observation": "sequencing_sample",
                "n_donors": None,
                "n_donors_status": "untestable_no_donor_key_on_disk",
                "fibrosis_manifest_value": manifest_fibrosis,
                "fibrosis_rederived_value": derived["fibrosis_stage"],
                "fibrosis_agrees": manifest_fibrosis == derived["fibrosis_stage"],
                "fibrosis_native_field": derived["fibrosis_native_field"],
                "fibrosis_native_value": derived["fibrosis_native_value"],
                "fibrosis_scale_native": FIBROSIS_SCALE_NATIVE[dataset],
                "cross_cohort_stage_comparable": CROSS_COHORT_STAGE_COMPARABLE[dataset],
                "nas_manifest_value": manifest_nas,
                "nas_rederived_value": derived["nas_score"],
                "nas_agrees": manifest_nas == derived["nas_score"],
                "nas_native_field": derived["nas_native_field"],
                "nas_native_value": derived["nas_native_value"],
                "sex_source": SEX_SOURCE[dataset],
                "source_table_role": (
                    "archived_at_risk" if dataset in AT_RISK_COHORTS else "pipeline_custom"
                ),
            }
        )

    manifest_runs = {(row["dataset"], row["sample_id"]) for row in manifest_rows}
    orphans: list[dict[str, Any]] = []
    for dataset, by_run in indexed.items():
        for run in sorted(by_run):
            if (dataset, run) not in manifest_runs:
                orphans.append(
                    {
                        "dataset": dataset,
                        "sra_run_accession": run,
                        "unit_of_observation": "sequencing_sample",
                        "present_in_source_table": True,
                        "present_in_manifest": False,
                        "yaml_excluded_samples": "empty_list",
                        "exclusion_reason": "untraced_downstream_of_metadata_harmonization",
                    }
                )
    return reproduction, orphans


def summarize(reproduction: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Per-cohort and total denominators, recomputed rather than copied from anywhere."""

    per_cohort: dict[str, dict[str, int]] = {}
    for row in reproduction:
        counts = per_cohort.setdefault(
            row["dataset"],
            {"n_samples": 0, "n_fibrosis_complete": 0, "n_nas_complete": 0,
             "n_both_complete": 0, "n_fibrosis_disagreements": 0, "n_nas_disagreements": 0},
        )
        counts["n_samples"] += 1
        has_fibrosis = row["fibrosis_rederived_value"] is not None
        has_nas = row["nas_rederived_value"] is not None
        counts["n_fibrosis_complete"] += int(has_fibrosis)
        counts["n_nas_complete"] += int(has_nas)
        counts["n_both_complete"] += int(has_fibrosis and has_nas)
        counts["n_fibrosis_disagreements"] += int(not row["fibrosis_agrees"])
        counts["n_nas_disagreements"] += int(not row["nas_agrees"])

    totals = {key: sum(c[key] for c in per_cohort.values())
              for key in next(iter(per_cohort.values()))}
    return {"per_cohort": per_cohort, "totals": totals}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-root", required=True,
                        help="RNA-seq/Human/Patient_Cohorts")
    parser.add_argument("--dataset-config", required=True,
                        help="config/human_datasets.yaml")
    parser.add_argument("--harmonize-script", required=True,
                        help="analysis/integration/scripts/00_harmonize_metadata.R")
    parser.add_argument("--manifest", required=True,
                        help="five_cohort_sample_manifest.tsv")
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output", required=True, help="staging directory")
    args = parser.parse_args(argv)

    cohort_root = Path(args.cohort_root).resolve(strict=True)
    output = Path(args.output)
    frozen_dir = output / "frozen"
    tables_dir = output / "tables"

    manifest_path = Path(args.manifest).resolve(strict=True)
    observed_manifest_sha = sha256_file(manifest_path)
    if observed_manifest_sha != args.manifest_sha256:
        raise ProvenanceError(
            "manifest digest changed since the freeze was planned: "
            f"expected {args.manifest_sha256}, observed {observed_manifest_sha}"
        )

    frozen: list[FrozenInput] = []
    source_tables: dict[str, list[dict[str, str]]] = {}
    for dataset, relative in COHORT_SOURCES.items():
        source = cohort_root / relative
        target = frozen_dir / "sra_run_tables" / f"{dataset}_SraRunTable.csv"
        frozen.append(freeze_input(source, target, role=f"sra_run_table_{dataset}"))
        source_tables[dataset] = read_csv_rows(target)

    frozen.append(freeze_input(Path(args.dataset_config).resolve(strict=True),
                               frozen_dir / "human_datasets.yaml",
                               role="dataset_registry"))
    frozen.append(freeze_input(Path(args.harmonize_script).resolve(strict=True),
                               frozen_dir / "00_harmonize_metadata.R",
                               role="harmonization_producer"))
    frozen.append(freeze_input(manifest_path,
                               frozen_dir / "five_cohort_sample_manifest.tsv",
                               role="substrate_manifest"))

    manifest_rows = read_tsv_rows(frozen_dir / "five_cohort_sample_manifest.tsv")
    reproduction, orphans = build_reproduction(manifest_rows, source_tables)
    summary = summarize(reproduction)

    write_tsv(tables_dir / "clinical_join_reproduction.tsv", REPRODUCTION_FIELDS, reproduction)
    write_tsv(
        tables_dir / "source_rows_absent_from_manifest.tsv",
        ("dataset", "sra_run_accession", "unit_of_observation", "present_in_source_table",
         "present_in_manifest", "yaml_excluded_samples", "exclusion_reason"),
        orphans,
    )
    write_tsv(
        tables_dir / "at_risk_inputs.tsv",
        ("dataset", "source_path", "frozen_path", "sha256", "include_in_mega",
         "deletion_consequence"),
        [
            {
                "dataset": dataset,
                "source_path": str(cohort_root / COHORT_SOURCES[dataset]),
                "frozen_path": f"frozen/sra_run_tables/{dataset}_SraRunTable.csv",
                "sha256": next(f.sha256 for f in frozen
                               if f.role == f"sra_run_table_{dataset}"),
                "include_in_mega": True,
                "deletion_consequence": (
                    "hard_stop_at_bg002_fail_closed_gate_00_harmonize_metadata_R_66_72"
                ),
            }
            for dataset in AT_RISK_COHORTS
        ],
    )
    write_tsv(
        output / "provenance" / "frozen_inputs.tsv",
        ("role", "source_path", "frozen_path", "sha256"),
        [
            {
                "role": item.role,
                "source_path": str(item.source_path),
                "frozen_path": str(item.frozen_path.relative_to(output)),
                "sha256": item.sha256,
            }
            for item in frozen
        ],
    )

    disagreements = summary["totals"]["n_fibrosis_disagreements"] + \
        summary["totals"]["n_nas_disagreements"]
    audit = {
        "unit_of_observation": "sequencing_sample",
        "donor_level_claim_made": False,
        "n_donors_status": "untestable_no_donor_key_on_disk",
        "join_key": "Run",
        "n_manifest_samples": len(manifest_rows),
        "n_reproduced": len(reproduction),
        "n_disagreements": disagreements,
        "n_source_rows_absent_from_manifest": len(orphans),
        "manifest_sha256": observed_manifest_sha,
        "summary": summary,
    }
    (output / "provenance").mkdir(parents=True, exist_ok=True)
    (output / "provenance" / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    if disagreements:
        raise ProvenanceError(
            f"{disagreements} clinical values could not be re-derived from the frozen "
            "source tables; the manifest join is not reproducible"
        )
    print(json.dumps({"n_reproduced": len(reproduction), "n_disagreements": 0,
                      "totals": summary["totals"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
