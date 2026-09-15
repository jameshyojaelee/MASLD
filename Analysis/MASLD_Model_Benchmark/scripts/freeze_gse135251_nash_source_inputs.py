#!/usr/bin/env python3
"""Freeze a hashed GSE135251 copy inside the benchmark as the NASH training source.

GSE135251 lives outside the benchmark tree and is read-only there.  This step
copies the exact bytes it read, records their SHA-256, and derives the source
fixture the transfer fit consumes.  Nothing is written back to the source tree
and GSE83452 is never opened here: this module has no argument that names it.

Three defects in that source tree are handled explicitly rather than absorbed.

1.  Every metadata file under that directory was written by R
    ``write.table(row.names = TRUE)`` without ``col.names = NA``.  The header
    carries N fields and every data row carries N + 1, because data field 0 is
    an unnamed rowname.  A naive ``csv.DictReader`` shifts every column by one
    and hands back the SRA run accession where ``Fibrosis_stage`` belongs.  The
    reader asserts the N + 1 invariant and fails closed, then cross-validates
    the parse against four independently declared census distributions.

2.  ``GSE135251_counts.tsv`` deposits 216 sample columns but the metadata
    describes 180.  The 36 columns with no metadata row have no outcome and are
    dropped, recorded, and never counted as participants.

3.  Its row key is a symbol/Ensembl mixture: 41,417 HGNC-style symbols and
    35,661 unversioned ENSG identifiers.  GSE83452 is keyed by Ensembl, so the
    join needs a crosswalk.  It is built from the GENCODE v49 GTF already
    hash-pinned in ``config/evaluation/microarray_transfer_preprocessing.toml``,
    kept to one-to-one mappings in both directions, and every unresolved row is
    retained as an explicit audit state rather than silently dropped.

The source NASH label is the deposited ``group_in_paper`` histology group.  Two
arms are emitted side by side and both are named before any GSE83452 outcome is
joined:

  ``nash_vs_nafl``            131 NASH against 41 NAFL.  The frozen label
                              requirements' ``prohibited_negative_substitutions``
                              bars ``healthy_control`` and ``normal_liver`` from
                              the negative class, and GSE83452's own negative
                              semantics say its no-NASH arm "is not a
                              healthy-control label", so the eight
                              ``disease = Control`` participants are excluded.

  ``nash_vs_nafl_plus_control`` 131 NASH against 49.  Recorded because it was
                              requested, not eligible under the criteria because it is the
                              substitution the requirements prohibit.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Any, Mapping, Sequence

import numpy as np


SERIES = "GSE135251"
EXPECTED_METADATA_ROWS = 180
EXPECTED_COUNTS_COLUMNS = 216
EXPECTED_COUNTS_ROWS = 77_078
GENCODE_RELEASE = "GENCODE_v49"

# Declared independently of this script, in the lane briefing.  Reproducing all
# four is what proves the N+1 parse landed on the right columns rather than one
# column over; `group_in_paper` and `Fibrosis_stage` also have to agree with each
# other, which a shifted parse cannot fake.
EXPECTED_CENSUS: dict[str, dict[str, int]] = {
    "group_in_paper": {
        "control": 8,
        "NAFL": 41,
        "NASH_F0-F1": 28,
        "NASH_F2": 47,
        "NASH_F3": 44,
        "NASH_F4": 12,
    },
    "Fibrosis_stage": {"0": 35, "1": 41, "2": 48, "3": 44, "4": 12},
    "nas_score": {
        "0": 8,
        "1": 7,
        "2": 16,
        "3": 22,
        "4": 34,
        "5": 42,
        "6": 32,
        "7": 13,
        "8": 6,
    },
    "disease": {"Control": 8, "NAFLD": 172},
}

NASH_GROUPS = ("NASH_F0-F1", "NASH_F2", "NASH_F3", "NASH_F4")
NAFL_GROUP = "NAFL"
CONTROL_GROUP = "control"

ARMS: dict[str, dict[str, Any]] = {
    "nash_vs_nafl": {
        "positive_groups": list(NASH_GROUPS),
        "negative_groups": [NAFL_GROUP],
        "excluded_groups": [CONTROL_GROUP],
        "expected_class_counts": {"no_nash": 41, "nash": 131},
        "gate_eligible_arm": True,
        "rationale": (
            "the frozen label contract bars healthy_control and normal_liver "
            "from the negative class, and GSE83452's deposited no-NASH arm is "
            "explicitly not a healthy-control label"
        ),
    },
    "nash_vs_nafl_plus_control": {
        "positive_groups": list(NASH_GROUPS),
        "negative_groups": [NAFL_GROUP, CONTROL_GROUP],
        "excluded_groups": [],
        "expected_class_counts": {"no_nash": 49, "nash": 131},
        "gate_eligible_arm": False,
        "rationale": (
            "recorded because it was requested; gate-ineligible because folding "
            "a healthy-control arm into no-NASH is the substitution "
            "prohibited_negative_substitutions names"
        ),
    },
}

AGE_FIELD_NAMES = frozenset(
    {"age", "age_years", "age_at_biopsy", "ageatbiopsy", "patient_age", "donor_age"}
)
SEX_FIELD_NAMES = frozenset({"sex", "gender", "patient_sex", "donor_sex"})

GENE_LINE = re.compile(r'\bgene_id "([^"]+)"')
GENE_NAME = re.compile(r'\bgene_name "([^"]+)"')


class SourceFreezeError(RuntimeError):
    """Raised when the GSE135251 freeze would misparse or misattribute a label."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# The N+1 metadata reader.
# ---------------------------------------------------------------------------


def read_rowname_offset_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read an R ``write.table(row.names = TRUE)`` table without shifting columns.

    Header field i names data field i + 1.  Data field 0 is the unnamed rowname.
    Any row that does not carry exactly N + 1 fields means the assumption behind
    this reader is wrong for this file, and that is an error, not a warning.
    """

    text = path.read_text(encoding="utf-8")
    lines = [line for line in text.split("\n") if line != ""]
    if len(lines) < 2:
        raise SourceFreezeError(f"metadata table has no data rows: {path}")
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for index, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header) + 1:
            raise SourceFreezeError(
                f"{path} line {index} carries {len(fields)} fields, not "
                f"{len(header) + 1}; the R rowname offset does not hold"
            )
        record = {header[i]: fields[i + 1] for i in range(len(header))}
        record["__rowname"] = fields[0]
        rows.append(record)
    return header, rows


def assert_declared_census(rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    audit: dict[str, Any] = {}
    for field, expected in EXPECTED_CENSUS.items():
        observed = Counter(row[field] for row in rows)
        if dict(observed) != expected:
            raise SourceFreezeError(
                f"{field} census differs from the declared distribution: "
                f"{dict(sorted(observed.items()))}"
            )
        audit[field] = {key: int(observed[key]) for key in sorted(observed)}
    # group_in_paper and Fibrosis_stage must agree with each other.  A parse that
    # is one column off cannot satisfy this and the four censuses at once.
    for row in rows:
        group = row["group_in_paper"]
        stage = row["Fibrosis_stage"]
        if group == "NASH_F2" and stage != "2":
            raise SourceFreezeError("NASH_F2 does not carry fibrosis stage 2")
        if group == "NASH_F3" and stage != "3":
            raise SourceFreezeError("NASH_F3 does not carry fibrosis stage 3")
        if group == "NASH_F4" and stage != "4":
            raise SourceFreezeError("NASH_F4 does not carry fibrosis stage 4")
        if group == "NASH_F0-F1" and stage not in {"0", "1"}:
            raise SourceFreezeError("NASH_F0-F1 does not carry fibrosis stage 0 or 1")
    audit["group_in_paper_agrees_with_fibrosis_stage"] = True
    # Matched exactly, never as a substring: "Fibrosis_stage" and "Stage" both
    # contain "age", and a substring test would report an age covariate that is
    # not there and silently reverse the age_sex_logistic unfittable finding.
    present = {field.strip().lower().replace(" ", "_") for field in rows[0]}
    audit["age_field_present"] = bool(present & AGE_FIELD_NAMES)
    audit["sex_field_present"] = bool(present & SEX_FIELD_NAMES)
    return audit


# ---------------------------------------------------------------------------
# The GENCODE v49 one-to-one crosswalk.
# ---------------------------------------------------------------------------


def read_gencode_genes(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Return gene_id -> gene_name and gene_name -> sorted gene_ids, unversioned."""

    observed: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = fields[8]
            gene_id = GENE_LINE.search(attributes)
            gene_name = GENE_NAME.search(attributes)
            if gene_id is None:
                continue
            stable = gene_id.group(1).split(".", 1)[0]
            observed[stable].add(gene_name.group(1) if gene_name is not None else "")
    # A stable ID that carries two different symbols across scaffolds is the
    # ambiguity the preprocessing requirements calls join_unresolved.  It keeps its
    # place in the ID census so an ENSG-keyed row still resolves, but it names no
    # symbol, so no symbol-keyed row can resolve onto it.
    id_to_name = {
        stable: (next(iter(names)) if len(names) == 1 else "")
        for stable, names in observed.items()
    }
    name_to_ids: dict[str, set[str]] = defaultdict(set)
    for stable, name in id_to_name.items():
        if name:
            name_to_ids[name].add(stable)
    return id_to_name, {name: sorted(ids) for name, ids in name_to_ids.items()}


def build_crosswalk(
    *, row_keys: Sequence[str], id_to_name: Mapping[str, str],
    name_to_ids: Mapping[str, Sequence[str]],
) -> tuple[dict[str, str], list[dict[str, str]]]:
    """Resolve each counts row key to at most one GENCODE v49 stable gene ID."""

    provisional: dict[str, str] = {}
    states: list[dict[str, str]] = []
    for key in row_keys:
        if key.startswith("ENSG"):
            if "." in key:
                states.append({"row_key": key, "mapping_state": "versioned_identifier"})
                continue
            if key in id_to_name:
                provisional[key] = key
                states.append(
                    {"row_key": key, "mapping_state": "ensembl_direct", "stable_gene_id": key}
                )
            else:
                states.append(
                    {"row_key": key, "mapping_state": "absent_from_gencode_v49"}
                )
            continue
        candidates = name_to_ids.get(key, ())
        if len(candidates) == 1:
            stable = candidates[0]
            if id_to_name.get(stable) != key:
                states.append({"row_key": key, "mapping_state": "join_unresolved"})
                continue
            provisional[key] = stable
            states.append(
                {"row_key": key, "mapping_state": "symbol_one_to_one", "stable_gene_id": stable}
            )
        elif len(candidates) > 1:
            states.append({"row_key": key, "mapping_state": "ambiguous_symbol"})
        else:
            states.append({"row_key": key, "mapping_state": "unmapped_symbol"})
    # Two row keys landing on one gene are a collision, not a gene.  Both go to
    # join_unresolved; neither is picked over the other.
    claims = Counter(provisional.values())
    resolved = {
        key: stable for key, stable in provisional.items() if claims[stable] == 1
    }
    collided = {key for key, stable in provisional.items() if claims[stable] > 1}
    for record in states:
        if record["row_key"] in collided:
            record["mapping_state"] = "join_unresolved_collision"
            record.pop("stable_gene_id", None)
    for record in states:
        record.setdefault("stable_gene_id", "")
    return resolved, states


# ---------------------------------------------------------------------------
# Fixture assembly.
# ---------------------------------------------------------------------------


def read_counts(
    path: Path, *, wanted_columns: Sequence[str], resolved: Mapping[str, str]
) -> tuple[list[str], np.ndarray, int]:
    """Read the counts matrix onto the resolved gene axis for the wanted columns."""

    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        columns = header[1:]
        if len(columns) != EXPECTED_COUNTS_COLUMNS:
            raise SourceFreezeError("counts matrix column census differs")
        position = {value: index for index, value in enumerate(columns)}
        if len(position) != len(columns):
            raise SourceFreezeError("counts matrix repeats a sample column")
        missing = [value for value in wanted_columns if value not in position]
        if missing:
            raise SourceFreezeError(f"counts matrix lacks metadata samples: {missing[:5]}")
        take = [position[value] for value in wanted_columns]
        collected: dict[str, list[float]] = {}
        rows_seen = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(columns) + 1:
                raise SourceFreezeError("counts matrix row is ragged")
            rows_seen += 1
            stable = resolved.get(fields[0])
            if stable is None:
                continue
            collected[stable] = [float(fields[1 + index]) for index in take]
    if rows_seen != EXPECTED_COUNTS_ROWS:
        raise SourceFreezeError("counts matrix row census differs")
    gene_ids = sorted(collected)
    matrix = np.asarray(
        [collected[gene_id] for gene_id in gene_ids], dtype=np.float64
    ).T
    if matrix.shape != (len(wanted_columns), len(gene_ids)):
        raise SourceFreezeError("assembled source matrix shape differs")
    if not np.all(np.isfinite(matrix)) or np.any(matrix < 0):
        raise SourceFreezeError("source counts are not finite and nonnegative")
    return gene_ids, matrix, rows_seen


def assign_outer_folds(labels: np.ndarray, *, folds: int, seed: int) -> np.ndarray:
    """Stratified outer folds over source participants, seeded and deterministic.

    GSE135251 deposits one run, one BioSample and one GEO sample per participant,
    so a participant cannot straddle a fold here; the assertion is in the caller.
    """

    rng = np.random.default_rng(seed)
    assignment = np.full(len(labels), -1, dtype=np.int64)
    for value in (0, 1):
        indices = np.flatnonzero(labels == value)
        shuffled = rng.permutation(indices)
        assignment[shuffled] = np.arange(len(shuffled)) % folds
    if np.any(assignment < 0):
        raise SourceFreezeError("outer-fold assignment is incomplete")
    for fold in range(folds):
        if set(int(v) for v in labels[assignment != fold]) != {0, 1}:
            raise SourceFreezeError("an outer-fold training partition lacks a class")
        if len(np.flatnonzero(assignment == fold)) == 0:
            raise SourceFreezeError("an outer fold is empty")
    return assignment


def run(
    *,
    counts_path: Path,
    metadata_path: Path,
    gtf_path: Path,
    gtf_sha256: str,
    output: Path,
    outer_folds: int,
    fold_seed: int,
) -> dict[str, Any]:
    from masld_bench.artifacts import freeze_tree

    if output.exists():
        raise SourceFreezeError(f"refusing to overwrite source freeze: {output}")
    counts_sha256 = sha256_file(counts_path)
    metadata_sha256 = sha256_file(metadata_path)
    observed_gtf_sha256 = sha256_file(gtf_path)
    if observed_gtf_sha256 != gtf_sha256:
        raise SourceFreezeError(
            "GENCODE v49 GTF SHA-256 differs from the frozen preprocessing contract"
        )

    header, metadata_rows = read_rowname_offset_tsv(metadata_path)
    if len(metadata_rows) != EXPECTED_METADATA_ROWS:
        raise SourceFreezeError("metadata row census differs")
    census = assert_declared_census(metadata_rows)
    if census["age_field_present"] or census["sex_field_present"]:
        raise SourceFreezeError(
            "an age or sex field appeared in the source metadata; the "
            "age_sex_logistic unfittable declaration must be revisited"
        )
    participant_ids = [row["GEO_Accession (exp)"] for row in metadata_rows]
    runs = [row["Run"] for row in metadata_rows]
    if len(set(participant_ids)) != len(participant_ids):
        raise SourceFreezeError("a source participant identifier repeats")
    if len(set(runs)) != len(runs):
        raise SourceFreezeError("a source run accession repeats")
    order = sorted(range(len(metadata_rows)), key=lambda index: participant_ids[index])
    metadata_rows = [metadata_rows[index] for index in order]
    participant_ids = [participant_ids[index] for index in order]
    runs = [runs[index] for index in order]

    id_to_name, name_to_ids = read_gencode_genes(gtf_path)
    with counts_path.open(encoding="utf-8", newline="") as handle:
        counts_header = handle.readline().rstrip("\n").split("\t")
        row_keys = [line.split("\t", 1)[0] for line in handle if line.strip()]
    if len(row_keys) != len(set(row_keys)):
        raise SourceFreezeError("counts matrix repeats a row key")
    resolved, mapping_states = build_crosswalk(
        row_keys=row_keys, id_to_name=id_to_name, name_to_ids=name_to_ids
    )
    gene_ids, matrix, rows_seen = read_counts(
        counts_path, wanted_columns=runs, resolved=resolved
    )

    output.mkdir(parents=True)
    raw_root = output / "raw"
    molecular = output / "molecular"
    outcomes = output / "outcomes"
    folds_root = output / "folds"
    for directory in (raw_root, molecular, outcomes, folds_root):
        directory.mkdir()

    shutil.copyfile(counts_path, raw_root / counts_path.name)
    shutil.copyfile(metadata_path, raw_root / metadata_path.name)
    for name, expected in (
        (counts_path.name, counts_sha256),
        (metadata_path.name, metadata_sha256),
    ):
        if sha256_file(raw_root / name) != expected:
            raise SourceFreezeError(f"frozen copy of {name} differs from the source")

    np.save(molecular / "rna_values.npy", matrix, allow_pickle=False)
    write_tsv(
        molecular / "rna_feature_axis.tsv",
        ("feature_index", "stable_gene_id"),
        [
            {"feature_index": index, "stable_gene_id": gene_id}
            for index, gene_id in enumerate(gene_ids)
        ],
    )
    write_tsv(
        molecular / "participant_axis.tsv",
        ("participant_id", "run_accession"),
        [
            {"participant_id": participant_id, "run_accession": run_accession}
            for participant_id, run_accession in zip(participant_ids, runs, strict=True)
        ],
    )
    write_tsv(
        molecular / "gene_mapping_states.tsv",
        ("row_key", "mapping_state", "stable_gene_id"),
        sorted(mapping_states, key=lambda record: record["row_key"]),
    )

    endpoint_rows: list[dict[str, Any]] = []
    for index, row in enumerate(metadata_rows):
        group = row["group_in_paper"]
        record: dict[str, Any] = {
            "participant_id": participant_ids[index],
            "group_in_paper": group,
            "disease": row["disease"],
            "fibrosis_stage": row["Fibrosis_stage"],
            "nas_score": row["nas_score"],
        }
        for arm_id, arm in ARMS.items():
            if group in arm["positive_groups"]:
                record[f"nash_status__{arm_id}"] = "nash"
            elif group in arm["negative_groups"]:
                record[f"nash_status__{arm_id}"] = "no_nash"
            else:
                record[f"nash_status__{arm_id}"] = "excluded"
        endpoint_rows.append(record)
    endpoint_fields = (
        "participant_id",
        "group_in_paper",
        "disease",
        "fibrosis_stage",
        "nas_score",
    ) + tuple(f"nash_status__{arm_id}" for arm_id in sorted(ARMS))
    write_tsv(outcomes / "participant_endpoints.tsv", endpoint_fields, endpoint_rows)

    arm_audit: dict[str, Any] = {}
    fold_rows: list[dict[str, Any]] = []
    fold_by_arm: dict[str, dict[str, int]] = {}
    for arm_id, arm in sorted(ARMS.items()):
        values = [row[f"nash_status__{arm_id}"] for row in endpoint_rows]
        counts = Counter(values)
        observed = {
            "no_nash": int(counts.get("no_nash", 0)),
            "nash": int(counts.get("nash", 0)),
        }
        if observed != arm["expected_class_counts"]:
            raise SourceFreezeError(
                f"{arm_id} class census differs: {observed}"
            )
        keep = [index for index, value in enumerate(values) if value != "excluded"]
        labels = np.asarray(
            [1 if values[index] == "nash" else 0 for index in keep], dtype=np.int64
        )
        assignment = assign_outer_folds(labels, folds=outer_folds, seed=fold_seed)
        fold_by_arm[arm_id] = {
            participant_ids[keep[position]]: int(assignment[position])
            for position in range(len(keep))
        }
        arm_audit[arm_id] = {
            "arm_id": arm_id,
            "positive_groups": arm["positive_groups"],
            "negative_groups": arm["negative_groups"],
            "excluded_groups": arm["excluded_groups"],
            "gate_eligible_arm": bool(arm["gate_eligible_arm"]),
            "rationale": arm["rationale"],
            "class_counts": observed,
            "fitting_participants": int(len(keep)),
            "excluded_participants": int(len(values) - len(keep)),
            "training_prevalence_nash": float(np.mean(labels)),
            "outer_folds": int(outer_folds),
            "outer_fold_seed": int(fold_seed),
            "outer_fold_sizes": {
                str(fold): int(np.sum(assignment == fold)) for fold in range(outer_folds)
            },
        }
    for participant_id in participant_ids:
        record: dict[str, Any] = {"participant_id": participant_id}
        for arm_id in sorted(ARMS):
            fold = fold_by_arm[arm_id].get(participant_id)
            record[f"outer_fold__{arm_id}"] = "excluded" if fold is None else str(fold)
        fold_rows.append(record)
    write_tsv(
        folds_root / "participant_outer_folds.tsv",
        ("participant_id",) + tuple(f"outer_fold__{arm_id}" for arm_id in sorted(ARMS)),
        fold_rows,
    )

    state_counts = Counter(record["mapping_state"] for record in mapping_states)
    receipt = {
        "schema_version": "masld-bench-gse135251-nash-source-freeze-v1",
        "status": "pass_source_frozen_hashed_copy",
        "source_series": SERIES,
        "external_series_read": False,
        "external_series": "GSE83452",
        "source_tree": str(counts_path.parent),
        "source_tree_is_read_only_to_this_step": True,
        "source_files_sha256": {
            counts_path.name: counts_sha256,
            metadata_path.name: metadata_sha256,
        },
        "gencode_gtf": str(gtf_path),
        "gencode_gtf_sha256": observed_gtf_sha256,
        "gencode_release": GENCODE_RELEASE,
        "metadata_rowname_offset_invariant": "header_fields_plus_one",
        "metadata_header_fields": len(header),
        "metadata_rows": len(metadata_rows),
        "declared_census_reproduced": census,
        "counts_columns_deposited": EXPECTED_COUNTS_COLUMNS,
        "counts_columns_without_metadata": EXPECTED_COUNTS_COLUMNS
        - len(metadata_rows),
        "counts_rows_deposited": rows_seen,
        "counts_row_key_states": {key: int(state_counts[key]) for key in sorted(state_counts)},
        "source_genes_resolved": len(gene_ids),
        "participants": len(participant_ids),
        "one_record_per_participant": True,
        "age_available_in_source": False,
        "sex_available_in_source": False,
        "age_sex_logistic_fittable_on_this_source": False,
        "age_sex_logistic_unfittable_reason": (
            "GSE135251 deposits no age and no recorded sex field; the TaskSpec "
            "forbids inferring a missing clinical field, so an age-sex-only lane "
            "cannot be fit on this training source"
        ),
        "arms": arm_audit,
        "gate_eligible_arm_id": "nash_vs_nafl",
        "source_roster_admission_state": "proposed_not_approved",
        "source_roster_admission_note": (
            "microarray_transfer_label_semantics.toml carries "
            "eligible_source_roster = [] and transfer_training_allowed_now = "
            "false; this freeze proposes GSE135251 as a training source and does "
            "not record it as admitted"
        ),
        "external_labels_read": False,
        "external_expression_values_read": False,
    }
    with (output / "source_freeze_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    for directory, artifact_class in (
        (raw_root, "gse135251_nash_source_raw_copy"),
        (molecular, "gse135251_nash_source_molecular"),
        (outcomes, "gse135251_nash_source_outcomes"),
        (folds_root, "gse135251_nash_source_folds"),
    ):
        freeze_tree(
            directory,
            {
                "artifact_class": artifact_class,
                "source_series": SERIES,
                "participants": len(participant_ids),
                "external_labels_read": False,
                "external_expression_values_read": False,
                "status": "passed",
            },
        )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--gencode-gtf", required=True, type=Path)
    parser.add_argument("--gencode-gtf-sha256", required=True)
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--fold-seed", type=int, default=1701)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run(
        counts_path=arguments.counts,
        metadata_path=arguments.metadata,
        gtf_path=arguments.gencode_gtf,
        gtf_sha256=arguments.gencode_gtf_sha256,
        output=arguments.output,
        outer_folds=arguments.outer_folds,
        fold_seed=arguments.fold_seed,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
