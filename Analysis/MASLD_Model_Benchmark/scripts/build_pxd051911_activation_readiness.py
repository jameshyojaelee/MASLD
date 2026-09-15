#!/usr/bin/env python3
"""Audit the PXD051911 activation requirements without inferring an RNA edge.

The PRIDE deposit ships its own ``checksum.txt`` manifest of SHA-1 digests.
That manifest, not a project-side record, is the provenance authority here:
every staged matrix is authenticated against it before any topology claim is
made.  Liver, scWAT, and oWAT stay assay-native and separable; no cross-tissue
value is imputed and no RNA pairing is created.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha1, sha256
import json
from pathlib import Path


# Digests transcribed from the deposit's own checksum.txt.  A mismatch means
# the local copy is not the published output file and must not be activated.
DEPOSIT_SHA1 = {
    "liver_protein_quant.txt": "810f33d2bbffc562fb9842df5ba75eff0bde0d15",
    "meta_data.txt": "f2b9940a67ff2d405bb2e643818bbd4fd2d5d8d5",
    "owat_protein_quant.txt": "11bef8905d77b36650f7e1b6cb5651329e814407",
    "scwat_protein_quant.txt": "802340721eac658538b4f000c2981992c4176735",
}

# Histology and participant fields the source supplies for every liver donor.
REQUIRED_COMPLETE_FIELDS = (
    "gender",
    "alder",
    "bmi",
    "kleiner_fibrosis_grade",
    "steatosis_score",
    "hepatocellular_ballooning_score",
    "lobular_inflammation_score",
    "nafld_activity_score",
    "saf_diagnosis",
)

# Present in the source but incomplete; it earns a typed mask, never a zero.
MASKED_FIELDS = ("fat_pct",)

TISSUES = {
    "liver": ("liver_protein_quant.txt", "liver_proteomics_filename"),
    "scwat": ("scwat_protein_quant.txt", "scWAT_proteomics_filename"),
    "owat": ("owat_protein_quant.txt", "oWAT_proteomics_filename"),
}

MISSING = {"NA", "", "NaN", "nan"}


class PXD051911ActivationError(RuntimeError):
    """Raised when the source differs from the audited activation requirements."""


def _digest(path: Path, algorithm: str) -> str:
    value = sha1() if algorithm == "sha1" else sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def observed(value: str | None) -> bool:
    return value is not None and value not in MISSING


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def quant_sample_columns(path: Path) -> list[str]:
    """Return the per-run columns, dropping the three annotation columns."""

    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header[:3] != ["ProteinAccessions", "Genes", "ProteinDescriptions"]:
        raise PXD051911ActivationError(f"{path.name} annotation columns differ")
    return header[3:]


def quant_row_count(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(1 for _ in handle) - 1


def tissue_join(
    rows: list[dict[str, str]], column: str, columns: list[str]
) -> dict[str, object]:
    """Join one tissue's quant columns to participants by source filename."""

    named = {r[column]: r for r in rows if observed(r.get(column))}
    quant = set(columns)
    matched = sorted(quant & set(named))
    return {
        "metadata_rows_with_filename": len(named),
        "quant_sample_columns": len(columns),
        "exact_filename_join": len(matched),
        "metadata_filenames_without_quant_column": sorted(set(named) - quant),
        "quant_columns_without_metadata_row": sorted(quant - set(named)),
        "participants": sorted({named[f]["patient_name"] for f in matched}),
    }


def build_audit(source: Path, metadata: Path) -> dict[str, object]:
    rows = read_metadata(metadata)
    joins: dict[str, dict[str, object]] = {}
    proteins: dict[str, int] = {}
    for tissue, (filename, column) in TISSUES.items():
        path = source / filename
        columns = quant_sample_columns(path)
        join = tissue_join(rows, column, columns)
        if join["metadata_filenames_without_quant_column"]:
            raise PXD051911ActivationError(
                f"{tissue} metadata filenames lack a quant column"
            )
        if join["quant_columns_without_metadata_row"]:
            raise PXD051911ActivationError(
                f"{tissue} quant columns lack an authoritative metadata row"
            )
        joins[tissue] = join
        proteins[tissue] = quant_row_count(path)

    liver = set(joins["liver"]["participants"])
    scwat = set(joins["scwat"]["participants"])
    owat = set(joins["owat"]["participants"])
    liver_rows = [
        r
        for r in rows
        if observed(r.get("liver_proteomics_filename"))
    ]

    # The trap in this cohort: meta_data.txt is visit-level, not person-level.
    # The plasma arm repeats people across visits, so a task built from the
    # apparently clean plasma matrix would pseudoreplicate.  The liver arm does
    # not repeat anyone.  Both facts are derived here, never assumed.
    plasma_rows = [
        r
        for r in rows
        if observed(r.get("plasma_proteomics_filename"))
        or observed(r.get("plasma_follow_up_proteomics_filename"))
    ]
    plasma_people = {r["patient_name"] for r in plasma_rows}
    donor_safety = {
        "metadata_grain": "one row per participant-visit, not per participant",
        "metadata_rows": len(rows),
        "metadata_distinct_participants": len({r["patient_name"] for r in rows}),
        "metadata_sample_groups": dict(
            sorted(Counter(r["sample_group"] for r in rows).items())
        ),
        "liver_arm": {
            "rows": len(liver_rows),
            "distinct_participants": len(liver),
            "rows_per_participant": "1:1",
            "repeat_visit_rows": len(liver_rows) - len(liver),
            "sample_groups": dict(
                sorted(Counter(r["sample_group"] for r in liver_rows).items())
            ),
            "donor_safe": len(liver_rows) == len(liver),
        },
        "plasma_arm": {
            "rows": len(plasma_rows),
            "distinct_participants": len(plasma_people),
            "repeat_visit_rows": len(plasma_rows) - len(plasma_people),
            "donor_safe": len(plasma_rows) == len(plasma_people),
            "warning": "The plasma arm is longitudinal. Using it without "
            "collapsing to person_key would treat repeat visits from one "
            "participant as independent donors. The activated liver arm is "
            "unaffected; this cohort is not safe tissue-by-tissue by default.",
        },
    }

    # One liver aliquot per participant is what makes a donor-safe split
    # possible at all; repeat visits would make these technical replicates.
    groups = Counter(r["sample_group"] for r in liver_rows)
    if set(groups) != {"initial_sample"}:
        raise PXD051911ActivationError("liver participants span repeat visits")
    if len(liver) != len(liver_rows):
        raise PXD051911ActivationError("liver participant key is not unique")

    labels: dict[str, object] = {}
    for field in REQUIRED_COMPLETE_FIELDS:
        values = [r[field] for r in liver_rows]
        if any(not observed(v) for v in values):
            raise PXD051911ActivationError(f"{field} is incomplete on liver donors")
        labels[field] = {
            "state": "observed",
            "missing": 0,
            "distribution": dict(sorted(Counter(values).items())),
        }
    for field in MASKED_FIELDS:
        values = [r[field] for r in liver_rows]
        labels[field] = {
            "state": "partially_observed_typed_mask",
            "missing": sum(1 for v in values if not observed(v)),
            "encoded_as_zero": False,
        }

    return {
        "donor_safety": donor_safety,
        "joins": joins,
        "proteins_quantified": proteins,
        "liver_participants": len(liver),
        "liver_sample_groups": dict(groups),
        "cross_tissue": {
            "liver_and_scwat": len(liver & scwat),
            "liver_and_owat": len(liver & owat),
            "liver_and_both_adipose": len(liver & scwat & owat),
            "scwat_only": len(scwat - liver),
            "owat_only": len(owat - liver),
        },
        "labels": labels,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--deposit-checksum", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise PXD051911ActivationError("activation readiness output exists")

    # Authenticate against the deposit's own manifest, then against the
    # transcribed constants, so neither alone can silently drift.
    published: dict[str, str] = {}
    for line in args.deposit_checksum.read_text(encoding="utf-8").splitlines():
        if "\t" not in line or line.startswith("#"):
            continue
        path_text, digest = line.rsplit("\t", 1)
        published[path_text.replace("\\", "/").rsplit("/", 1)[-1]] = digest.strip()

    provenance: dict[str, object] = {}
    for filename, expected in DEPOSIT_SHA1.items():
        if published.get(filename) != expected:
            raise PXD051911ActivationError(
                f"deposit checksum.txt disagrees with the audited digest for {filename}"
            )
        path = args.source / filename
        actual = _digest(path, "sha1")
        if actual != expected:
            raise PXD051911ActivationError(
                f"{filename} SHA-1 {actual} differs from the published {expected}"
            )
        provenance[filename] = {
            "sha1_published": expected,
            "sha1_local": actual,
            "sha256": _digest(path, "sha256"),
            "size_bytes": path.stat().st_size,
            "authenticated_against": "deposit_checksum_txt",
        }

    audit = build_audit(args.source, args.source / "meta_data.txt")
    args.output.mkdir(parents=True, exist_ok=False)

    manifest = {
        "schema_version": "masld-bench-pxd051911-dataset-manifest-v1",
        "dataset_id": "pxd051911",
        "accession": "PXD051911",
        "biological_unit": "participant",
        "expected_biological_units": 58,
        # Stated at the top of the manifest on purpose. The metadata is
        # visit-level, so this cohort is not uniformly donor-safe: the liver
        # arm is 1:1 with participants and the plasma arm is not.
        "donor_safety": audit["donor_safety"],
        "activated_arm": {
            "tissues": ["liver_proteomics", "scwat_proteomics", "owat_proteomics"],
            "excluded_tissue": "plasma_proteomics",
            "exclusion_reason": "The plasma arm repeats participants across "
            "visits and would pseudoreplicate. It is not part of this "
            "activation and no plasma endpoint is registered.",
        },
        "rights": {
            "license": "Creative Commons Public Domain (CC0)",
            "license_source": "PRIDE project record",
            "access_tier": "public",
            "credential_gated": False,
            "redistribution_permitted": True,
        },
        "provenance": provenance,
        "source_authority": "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD051911",
        "reference": {
            "native_genome_build": "not_applicable",
            "search_database": "HomoSapiens_Oct19_Swiss_withiRT.fasta",
            "feature_axis": "uniprot_protein_accession",
            "coordinate_check": "not_applicable_proteomics",
        },
    }
    (args.output / "dataset_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    rows = read_metadata(args.source / "meta_data.txt")
    by_file = {
        r["liver_proteomics_filename"]: r
        for r in rows
        if observed(r.get("liver_proteomics_filename"))
    }
    scwat_by_patient = {
        r["patient_name"]: r["scWAT_proteomics_filename"]
        for r in rows
        if observed(r.get("scWAT_proteomics_filename"))
    }
    owat_by_patient = {
        r["patient_name"]: r["oWAT_proteomics_filename"]
        for r in rows
        if observed(r.get("oWAT_proteomics_filename"))
    }
    fieldnames = [
        "person_key",
        "source_patient_name",
        "liver_run",
        "scwat_run",
        "scwat_state",
        "owat_run",
        "owat_state",
        "sex",
        "age_years",
        "bmi",
        "fat_pct",
        "fat_pct_state",
        "kleiner_fibrosis_grade",
        "steatosis_score",
        "hepatocellular_ballooning_score",
        "lobular_inflammation_score",
        "nafld_activity_score",
        "saf_diagnosis",
    ]
    with (args.output / "participant_join.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for filename in sorted(by_file):
            row = by_file[filename]
            patient = row["patient_name"]
            scwat = scwat_by_patient.get(patient)
            owat = owat_by_patient.get(patient)
            writer.writerow(
                {
                    "person_key": f"pxd051911:{patient}",
                    "source_patient_name": patient,
                    "liver_run": filename,
                    "scwat_run": scwat or "",
                    "scwat_state": "observed" if scwat else "structurally_missing",
                    "owat_run": owat or "",
                    "owat_state": "observed" if owat else "structurally_missing",
                    "sex": row["gender"],
                    "age_years": row["alder"],
                    "bmi": row["bmi"],
                    "fat_pct": row["fat_pct"] if observed(row["fat_pct"]) else "",
                    "fat_pct_state": (
                        "observed" if observed(row["fat_pct"]) else "structurally_missing"
                    ),
                    "kleiner_fibrosis_grade": row["kleiner_fibrosis_grade"],
                    "steatosis_score": row["steatosis_score"],
                    "hepatocellular_ballooning_score": row[
                        "hepatocellular_ballooning_score"
                    ],
                    "lobular_inflammation_score": row["lobular_inflammation_score"],
                    "nafld_activity_score": row["nafld_activity_score"],
                    "saf_diagnosis": row["saf_diagnosis"],
                }
            )

    receipt = {
        "schema_version": "masld-bench-pxd051911-activation-readiness-v1",
        "dataset_id": "pxd051911",
        "status": "activation_ready_development_only",
        "unit_of_replication": "participant",
        "sealed_outcomes_read": False,
        "models_fit": [],
        "models_scored": [],
        "rna_protein_pair_inferred": False,
        "missing_encoded_as_zero": False,
        "split": {
            "group_key": "cohort_family_id+person_key",
            "donor_safe": True,
            "repeat_visit_rows_in_liver_arm": 0,
            "all_tissues_from_one_person_share_outer_role": True,
        },
        "champion_claim_allowed": False,
        "champion_blocker": (
            "protein_transport requires an external protein cohort and zero are "
            "eligible; PXD051911 alone is development-only."
        ),
        **audit,
    }
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
