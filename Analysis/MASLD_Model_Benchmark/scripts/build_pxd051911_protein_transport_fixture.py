"""Build the PXD051911 participant x protein fixture for the protein_transport task.

The fixture is source-native. Intensities are written exactly as deposited, with no
normalisation, no log transform, and no imputation, so every learned transform stays
inside the outer training participants where the TaskSpec requires it.

Missingness carries a typed four-code state and is never encoded as zero. A protein
absent from a tissue panel and a protein present in the panel but not quantified in a
participant are different facts with different error models; collapsing them would be
a scientific error, not a formatting one.
"""

from __future__ import annotations

import argparse
import csv
import json
from hashlib import sha256
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

FOLD_SEED = 20260821
OUTER_FOLDS = 5
COHORT_FAMILY_ID = "pxd051911"

# Typed protein state codes. Codes 1 and 2 are both `structurally_missing` in the
# MISSING_STATES vocabulary; they are kept apart because they arise for different
# reasons and a downstream model may legitimately treat them differently.
STATE_OBSERVED = 0
STATE_PROTEIN_NOT_IN_PANEL = 1
STATE_TISSUE_NOT_COLLECTED = 2
STATE_BELOW_QC = 3

STATE_CODES: Mapping[str, int] = {
    "observed": STATE_OBSERVED,
    "structurally_missing_protein_not_in_tissue_panel": STATE_PROTEIN_NOT_IN_PANEL,
    "structurally_missing_tissue_not_collected_for_participant": STATE_TISSUE_NOT_COLLECTED,
    "below_qc_not_quantified_in_this_participant": STATE_BELOW_QC,
}
STATE_TO_MISSING_STATE: Mapping[str, str] = {
    "observed": "observed",
    "structurally_missing_protein_not_in_tissue_panel": "structurally_missing",
    "structurally_missing_tissue_not_collected_for_participant": "structurally_missing",
    "below_qc_not_quantified_in_this_participant": "below_qc",
}

TISSUES = ("liver", "scwat", "owat")
QUANT_FILES = {
    "liver": "liver_protein_quant.txt",
    "scwat": "scwat_protein_quant.txt",
    "owat": "owat_protein_quant.txt",
}
RUN_COLUMN = {"liver": "liver_run", "scwat": "scwat_run", "owat": "owat_run"}

FIBROSIS_CODE = {"F0": 0, "F1": 1, "F2": 2, "F3": 3, "F4": 4}


class FixtureError(RuntimeError):
    """Raised when the fixture cannot be built exactly as contracted."""


def fold_index(unit_id: str, *, seed: int, outer_folds: int) -> int:
    """Mirror masld_bench.evaluators.rna_atac_development.fold_index exactly."""
    digest = sha256(f"{seed}\0{unit_id}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % outer_folds


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise FixtureError(f"refusing to write an empty table: {path}")
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_quant(path: Path) -> tuple[list[str], dict[str, str], dict[str, dict[str, float]]]:
    """Return (accessions, accession->gene symbol, run -> accession -> intensity).

    Missing cells are simply absent from the returned mapping. `NA` is the deposit's
    own missing token; an explicit numeric zero would be a real measurement and is
    rejected, because the deposit contains none and one appearing later would mean the
    upstream file changed meaning.
    """
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        runs = header[3:]
        accessions: list[str] = []
        genes: dict[str, str] = {}
        values: dict[str, dict[str, float]] = {run: {} for run in runs}
        for row in reader:
            accession, gene = row[0], row[1]
            accessions.append(accession)
            genes[accession] = gene
            for run, cell in zip(runs, row[3:]):
                if cell in ("NA", ""):
                    continue
                number = float(cell)
                if number == 0.0:
                    raise FixtureError(
                        f"explicit zero intensity in {path.name} at {accession}/{run}; "
                        "the deposit encodes missingness as NA only"
                    )
                values[run][accession] = number
    if len(accessions) != len(set(accessions)):
        raise FixtureError(f"duplicate protein accession in {path.name}")
    return accessions, genes, values


def build(source: Path, join_path: Path, output: Path) -> dict:
    join = read_tsv(join_path)
    if len(join) != 58:
        raise FixtureError(f"expected the 58-participant liver arm, found {len(join)}")

    participants = sorted(row["person_key"] for row in join)
    by_key = {row["person_key"]: row for row in join}
    index_of = {key: i for i, key in enumerate(participants)}
    n = len(participants)

    panels: dict[str, list[str]] = {}
    quants: dict[str, dict[str, dict[str, float]]] = {}
    genes: dict[str, str] = {}
    for tissue in TISSUES:
        accessions, tissue_genes, values = load_quant(source / QUANT_FILES[tissue])
        panels[tissue] = accessions
        quants[tissue] = values
        for accession, gene in tissue_genes.items():
            genes.setdefault(accession, gene)

    feature_axis = sorted(set().union(*(set(panels[t]) for t in TISSUES)))
    feature_of = {accession: i for i, accession in enumerate(feature_axis)}
    n_features = len(feature_axis)

    output.mkdir(parents=True, exist_ok=False)
    molecular = output / "molecular"
    folds_dir = output / "folds"
    molecular.mkdir()
    folds_dir.mkdir()

    write_tsv(
        molecular / "protein_feature_axis.tsv",
        [
            {
                "protein_feature_index": feature_of[accession],
                "protein_accession": accession,
                "source_gene_symbol": genes.get(accession, ""),
                "source_gene_symbol_state": (
                    "observed" if genes.get(accession, "") not in ("", "NA") else "structurally_missing"
                ),
                "in_liver_panel": str(accession in set(panels["liver"])).lower(),
                "in_scwat_panel": str(accession in set(panels["scwat"])).lower(),
                "in_owat_panel": str(accession in set(panels["owat"])).lower(),
            }
            for accession in feature_axis
        ],
    )

    receipt_tissues: dict[str, dict] = {}
    for tissue in TISSUES:
        panel = set(panels[tissue])
        values = np.full((n, n_features), np.nan, dtype=np.float64)
        states = np.full((n, n_features), STATE_BELOW_QC, dtype=np.uint8)

        not_in_panel = np.array(
            [accession not in panel for accession in feature_axis], dtype=np.bool_
        )
        states[:, not_in_panel] = STATE_PROTEIN_NOT_IN_PANEL

        collected = np.zeros(n, dtype=np.bool_)
        for key in participants:
            run = by_key[key][RUN_COLUMN[tissue]].strip()
            row = index_of[key]
            if not run:
                states[row, :] = STATE_TISSUE_NOT_COLLECTED
                continue
            if run not in quants[tissue]:
                raise FixtureError(f"{tissue} run {run!r} for {key} is absent from the quant matrix")
            collected[row] = True
            for accession, intensity in quants[tissue][run].items():
                column = feature_of[accession]
                values[row, column] = intensity
                states[row, column] = STATE_OBSERVED

        observed = states == STATE_OBSERVED
        if np.any(values[observed] == 0.0):
            raise FixtureError(f"{tissue}: an observed cell holds zero")
        if np.any(np.isnan(values[observed])):
            raise FixtureError(f"{tissue}: an observed cell holds NaN")
        if not np.all(np.isnan(values[~observed])):
            raise FixtureError(f"{tissue}: an unobserved cell holds a number")

        np.save(molecular / f"{tissue}_values.npy", values, allow_pickle=False)
        np.save(molecular / f"{tissue}_state.npy", states, allow_pickle=False)
        np.save(molecular / f"{tissue}_observed_mask.npy", observed, allow_pickle=False)
        np.save(molecular / f"{tissue}_participant_collected_mask.npy", collected, allow_pickle=False)

        receipt_tissues[tissue] = {
            "participants_with_tissue_collected": int(collected.sum()),
            "panel_proteins": len(panel),
            "observed_cells": int(observed.sum()),
            "below_qc_cells": int((states == STATE_BELOW_QC).sum()),
            "structurally_missing_protein_not_in_panel_cells": int(
                (states == STATE_PROTEIN_NOT_IN_PANEL).sum()
            ),
            "structurally_missing_tissue_not_collected_cells": int(
                (states == STATE_TISSUE_NOT_COLLECTED).sum()
            ),
            "within_panel_detection_rate": round(
                float(observed.sum()) / float(max(int(collected.sum()) * len(panel), 1)), 6
            ),
        }

    assignment = {key: fold_index(key, seed=FOLD_SEED, outer_folds=OUTER_FOLDS) for key in participants}
    write_tsv(
        folds_dir / "participant_outer_folds.tsv",
        [
            {
                "person_key": key,
                "cohort_family_id": COHORT_FAMILY_ID,
                "outer_fold": assignment[key],
            }
            for key in participants
        ],
    )
    fold_sizes = [sum(1 for k in participants if assignment[k] == f) for f in range(OUTER_FOLDS)]
    if sum(fold_sizes) != n:
        raise FixtureError("fold assignment lost a participant")

    endpoint_rows = []
    for key in participants:
        row = by_key[key]
        fat = row["fat_pct"].strip()
        endpoint_rows.append(
            {
                "person_key": key,
                "cohort_family_id": COHORT_FAMILY_ID,
                "outer_fold": assignment[key],
                "nafld_activity_score": int(row["nafld_activity_score"]),
                "steatosis_score": int(row["steatosis_score"]),
                "hepatocellular_ballooning_score": int(row["hepatocellular_ballooning_score"]),
                "lobular_inflammation_score": int(row["lobular_inflammation_score"]),
                "kleiner_fibrosis_grade": row["kleiner_fibrosis_grade"],
                "kleiner_fibrosis_code": FIBROSIS_CODE[row["kleiner_fibrosis_grade"]],
                "fibrosis_binary_f2_or_greater": int(FIBROSIS_CODE[row["kleiner_fibrosis_grade"]] >= 2),
                "saf_diagnosis": row["saf_diagnosis"],
                "sex": row["sex"],
                "age_years": int(row["age_years"]),
                "bmi": float(row["bmi"]),
                "fat_pct": fat,
                "fat_pct_state": "observed" if fat else "structurally_missing",
                "scwat_observation_state": row["scwat_state"],
                "owat_observation_state": row["owat_state"],
            }
        )
    write_tsv(output / "endpoints.tsv", endpoint_rows)

    for row in endpoint_rows:
        components = (
            row["steatosis_score"]
            + row["hepatocellular_ballooning_score"]
            + row["lobular_inflammation_score"]
        )
        if components != row["nafld_activity_score"]:
            raise FixtureError(f"NAS is not the component sum for {row['person_key']}")

    receipt = {
        "schema_version": "masld-bench-pxd051911-protein-transport-fixture-v1",
        "dataset_id": "pxd051911",
        "task_id": "protein_transport",
        "unit_of_replication": "participant",
        "participants": n,
        "participant_axis_is_liver_anchored": True,
        "participants_excluded_for_no_liver_anchor": sorted(
            {
                run.split("_")[-1].replace(".raw", "")
                for tissue in ("scwat", "owat")
                for run in quants[tissue]
                if run not in {by_key[k][RUN_COLUMN[tissue]].strip() for k in participants}
            }
        ),
        "feature_axis": {
            "union_proteins": n_features,
            "liver_panel": len(panels["liver"]),
            "scwat_panel": len(panels["scwat"]),
            "owat_panel": len(panels["owat"]),
            "feature_identity": "opaque_source_protein_accession",
            "gene_symbol_is_annotation_only_never_a_join_key": True,
        },
        "tissues": receipt_tissues,
        "values": {
            "source_native_raw_intensity": True,
            "normalisation_applied": False,
            "log_transform_applied": False,
            "imputation_applied": False,
            "note": "All normalisation, transformation and imputation must be fit inside the outer training participants by the consuming model.",
        },
        "missingness": {
            "state_codes": dict(STATE_CODES),
            "state_to_missing_state_vocabulary": dict(STATE_TO_MISSING_STATE),
            "missing_encoded_as_zero": False,
            "explicit_typed_state_matrix": True,
            "explicit_boolean_masks": True,
            "fat_pct_observed": sum(1 for r in endpoint_rows if r["fat_pct_state"] == "observed"),
            "fat_pct_structurally_missing": sum(
                1 for r in endpoint_rows if r["fat_pct_state"] == "structurally_missing"
            ),
        },
        "split": {
            "split_id": "pxd051911_participant_outer_v1",
            "group_key": "cohort_family_id+person_key",
            "outer_folds": OUTER_FOLDS,
            "seed": FOLD_SEED,
            "stratified": False,
            "fold_sizes": fold_sizes,
            "all_tissues_from_one_participant_share_one_fold": True,
        },
        "cross_tissue_scopes": {
            "S-SCWAT": int(sum(1 for k in participants if by_key[k]["scwat_state"] == "observed")),
            "S-OWAT": int(sum(1 for k in participants if by_key[k]["owat_state"] == "observed")),
            "pairing_level": "same_donor_different_tissue",
            "never_same_sample": True,
        },
        "endpoint_integrity": {
            "nas_equals_component_sum": True,
            "checked_participants": n,
        },
        "rna_protein_pair_inferred": False,
        "plasma_arm_used": False,
        "models_fit": [],
        "models_scored": [],
        "sealed_outcomes_read": False,
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--join", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipt = build(args.source, args.join, args.output)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
