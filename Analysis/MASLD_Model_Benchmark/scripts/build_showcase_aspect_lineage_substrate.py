"""Join, hash and freeze the aspect- and lineage-resolved showcase substrate.

Three arms, prepared and frozen only.  Nothing here fits a model or computes a
metric.

Arm A  GSE267145 aspect training substrate.  Bulk RNA joined to the deposited
       NASH-CRN component grades on the participant axis.
Arm B  GSE135251 external-check substrate.  The frozen n=180 source reused as
       deposited; the cohort deposits the NAS sum, not its components.
Arm C  GSE202379 aspect-and-lineage cross-check.  Donor-level SAF triples parsed
       from the primary GEO series matrix, joined to single-nucleus lineage
       counts through a primary SRA run table.

Load-bearing parsing rules, each one paid for already:

*   ``Analysis/Deconvolution/bulk`` metadata is R ``write.table(row.names=TRUE)``
    output: header N fields, data rows N + 1.  ``Analysis/Deconvolution/metadata``
    SRA tables are flat.  Both invariants are asserted per file, never assumed
    from a sibling.
*   BayesPrism proportions carry a named lineage header and an unnamed key
    column, so data rows hold one more field than the header.
*   The GSE202379 ``saf score`` field is never absent.  Ungraded means the value
    fails ``^S\\d+A\\d+F\\d+$``, never that the field is missing.
*   ``!Sample_characteristics_ch1`` rows are positionally ragged.  Every cell is
    split on its own ``key: value`` prefix; a row is never keyed by its first
    cell.
*   A row count is not a unit count.  Every census reports distinct units beside
    rows and names which one the arithmetic uses.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

SAF_TRIPLE = re.compile(r"^S(\d+)A(\d+)F(\d+)$")
ENSEMBL_UNVERSIONED = re.compile(r"^ENSG\d+$")
RUN_ACCESSION = re.compile(r"^SRR\d+$")
EXPERIMENT_ACCESSION = re.compile(r"^SRX\d+$")
SAMPLE_ACCESSION = re.compile(r"^GSM\d+$")

# The three NASH-CRN components that sum to the NAS.  Lobular necrosis is
# deposited beside them and is NOT one of them.
NAS_COMPONENTS = ("steatosis", "ballooning", "lobular_inflammation")
NOT_A_NAS_COMPONENT = "lobular_necrosis"
FOURTH_AXIS = "fibrosis"

# The REF_HUMAN BayesPrism lineage axis, in deposited header order.
REF_HUMAN_LINEAGES = (
    "Endothelial cells",
    "Hepatocytes",
    "Plasma cells",
    "T cells",
    "Cholangiocytes",
    "Fibroblasts",
    "Macrophages",
    "Circulating NK/NKT",
    "Resident NK",
    "Mono+mono derived cells",
    "Basophils",
    "B cells",
    "cDC1s",
    "cDC2s",
    "pDCs",
    "Neutrophils",
)


class SubstrateError(RuntimeError):
    """Raised when a joined substrate cannot be built exactly as specified."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            if len(row) != len(header):
                raise SubstrateError(f"{path}: row width {len(row)} != {len(header)}")
            handle.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def read_flat_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a table whose data rows carry exactly as many fields as the header.

    The sibling ``bulk/`` tables carry an R rowname offset and this one does not.
    Asserting which is which per file is the whole point: a reader that guesses
    returns the neighbouring column's value and never raises.
    """

    lines = [line for line in path.read_text(encoding="utf-8").split("\n") if line != ""]
    if len(lines) < 2:
        raise SubstrateError(f"table has no data rows: {path}")
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for index, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header):
            raise SubstrateError(
                f"{path} line {index} carries {len(fields)} fields, not "
                f"{len(header)}; this table is not flat"
            )
        rows.append(dict(zip(header, fields)))
    return header, rows


def read_bayesprism_proportions(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Read a BayesPrism proportions table.

    The header names the lineages only.  Data rows carry one extra leading field
    holding the unnamed sample key, so ``len(row) == len(header) + 1``.  A manual
    zip of header to row shifts every lineage by one and silently succeeds.
    """

    lines = [line for line in path.read_text(encoding="utf-8").split("\n") if line != ""]
    if len(lines) < 2:
        raise SubstrateError(f"proportions table has no data rows: {path}")
    lineages = lines[0].split("\t")
    keys: list[str] = []
    values: list[list[float]] = []
    for index, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(lineages) + 1:
            raise SubstrateError(
                f"{path} line {index} carries {len(fields)} fields, not "
                f"{len(lineages) + 1}; the unnamed key column is not present"
            )
        keys.append(fields[0])
        values.append([float(v) for v in fields[1:]])
    matrix = np.asarray(values, dtype=np.float64)
    if len(set(keys)) != len(keys):
        raise SubstrateError(f"{path}: sample keys are not unique")
    return lineages, keys, matrix


def parse_series_matrix_characteristics(
    path: Path,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    """Parse a GEO series matrix into per-GSM ``key: value`` characteristics.

    ``!Sample_characteristics_ch1`` rows are positionally ragged: a single row
    can carry ``gender`` in some columns and ``liver lobe`` in others.  Every
    cell is therefore split on its own prefix.  Keying a row by its first cell
    mixes fields and drops real values.
    """

    accessions: list[str] | None = None
    characteristic_rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                accessions = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristic_rows.append(
                    next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
                )
    if not accessions:
        raise SubstrateError(f"{path}: no !Sample_geo_accession row")
    if len(set(accessions)) != len(accessions):
        raise SubstrateError(f"{path}: sample accessions are not unique")
    fields: dict[str, dict[str, str]] = {acc: {} for acc in accessions}
    for row in characteristic_rows:
        if len(row) != len(accessions):
            raise SubstrateError(
                f"{path}: a characteristics row is {len(row)} wide, not "
                f"{len(accessions)}"
            )
        for accession, cell in zip(accessions, row):
            cell = cell.strip()
            if not cell or ":" not in cell:
                continue
            key, value = cell.split(":", 1)
            key = key.strip().lower().replace(" ", "_")
            value = value.strip()
            existing = fields[accession].get(key)
            if existing is not None and existing != value:
                raise SubstrateError(
                    f"{path}: {accession} carries conflicting {key!r} values"
                )
            fields[accession][key] = value
    return accessions, fields


def read_sra_run_info(path: Path) -> list[dict[str, str]]:
    """Read a headerless NCBI ``SraRunInfo.csv`` by column CONTENT, not by name.

    The file ships without a header row and quotes fields that contain embedded
    newlines, so a line-oriented split miscounts records.  Columns are located by
    the accession pattern every one of their cells matches.
    """

    with path.open(newline="", encoding="utf-8") as handle:
        records = list(csv.reader(handle))
    if not records:
        raise SubstrateError(f"{path}: empty run table")
    widths = {len(record) for record in records}
    if len(widths) != 1:
        raise SubstrateError(f"{path}: ragged run table, widths {sorted(widths)}")
    width = widths.pop()

    def columns_matching(pattern: re.Pattern[str]) -> list[int]:
        return [
            index
            for index in range(width)
            if all(pattern.fullmatch(record[index]) for record in records)
        ]

    run_columns = columns_matching(RUN_ACCESSION)
    experiment_columns = columns_matching(EXPERIMENT_ACCESSION)
    sample_columns = columns_matching(SAMPLE_ACCESSION)
    if len(run_columns) != 1 or not experiment_columns or not sample_columns:
        raise SubstrateError(
            f"{path}: could not locate exactly one run column and at least one "
            f"experiment and sample column (runs={run_columns}, "
            f"experiments={experiment_columns}, samples={sample_columns})"
        )
    for record in records:
        if len({record[index] for index in sample_columns}) != 1:
            raise SubstrateError(f"{path}: GSM-shaped columns disagree within a row")
    run, experiment, sample = (
        run_columns[0],
        experiment_columns[0],
        sample_columns[0],
    )
    rows = [
        {
            "run_accession": record[run],
            "experiment_accession": record[experiment],
            "sample_accession": record[sample],
        }
        for record in records
    ]
    if len({row["run_accession"] for row in rows}) != len(rows):
        raise SubstrateError(f"{path}: run accessions are not unique")
    return rows


def read_h5ad_obs_columns(path: Path, columns: Sequence[str]) -> dict[str, np.ndarray]:
    """Read named ``obs`` columns from an h5ad without materialising ``X``."""

    import h5py

    def decode(values: Iterable[Any]) -> np.ndarray:
        return np.asarray(
            [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values],
            dtype=object,
        )

    out: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        for column in columns:
            if column not in obs:
                raise SubstrateError(f"{path}: obs has no column {column!r}")
            node = obs[column]
            if isinstance(node, h5py.Group):
                if "categories" not in node or "codes" not in node:
                    raise SubstrateError(f"{path}: obs/{column} is not categorical")
                categories = decode(node["categories"][:])
                codes = node["codes"][:]
                if codes.min(initial=0) < 0:
                    raise SubstrateError(f"{path}: obs/{column} has unset codes")
                out[column] = categories[codes]
            else:
                out[column] = decode(node[:])
    lengths = {len(v) for v in out.values()}
    if len(lengths) != 1:
        raise SubstrateError(f"{path}: obs columns disagree in length")
    return out


def build_shared_gene_axis(
    arm_a_axis: Sequence[str], arm_b_axis: Sequence[str]
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Intersect two stable-gene-id axes.

    The intersection uses gene IDENTITY only.  No expression value and no
    outcome is read, so this filter cannot leak across an outer fold.  Any
    variance or abundance filter belongs inside outer training instead, which is
    why none is applied here.
    """

    for name, axis in (("arm A", arm_a_axis), ("arm B", arm_b_axis)):
        if len(set(axis)) != len(axis):
            raise SubstrateError(f"{name} gene axis carries duplicate stable ids")
        offenders = [g for g in axis if not ENSEMBL_UNVERSIONED.fullmatch(g)]
        if offenders:
            raise SubstrateError(
                f"{name} gene axis carries {len(offenders)} non-unversioned-ENSG "
                f"ids, first {offenders[:3]}"
            )
    a_index = {gene: i for i, gene in enumerate(arm_a_axis)}
    b_index = {gene: i for i, gene in enumerate(arm_b_axis)}
    shared = sorted(set(a_index) & set(b_index))
    if not shared:
        raise SubstrateError("arm A and arm B share no genes")
    return (
        shared,
        np.asarray([a_index[g] for g in shared], dtype=np.int64),
        np.asarray([b_index[g] for g in shared], dtype=np.int64),
    )


def build_arm_a(
    *,
    fixture: Path,
    activation: Path,
    output: Path,
    shared_genes: Sequence[str],
    shared_columns: np.ndarray,
) -> dict[str, Any]:
    molecular = fixture / "molecular"
    participant_header, participant_rows = read_flat_tsv(
        molecular / "participant_axis.tsv"
    )
    _, feature_rows = read_flat_tsv(molecular / "rna_feature_axis.tsv")
    _, fold_rows = read_flat_tsv(fixture / "folds" / "participant_outer_folds.tsv")
    _, endpoint_rows = read_flat_tsv(activation / "participant_endpoints.tsv")

    values = np.load(molecular / "rna_values.npy")
    if values.shape != (len(participant_rows), len(feature_rows)):
        raise SubstrateError(
            f"arm A matrix {values.shape} does not match "
            f"{len(participant_rows)} participants x {len(feature_rows)} features"
        )
    indices = [int(row["participant_index"]) for row in participant_rows]
    if indices != list(range(len(participant_rows))):
        raise SubstrateError("arm A participant_index is not 0..n-1 in row order")
    feature_indices = [int(row["rna_feature_index"]) for row in feature_rows]
    if feature_indices != list(range(len(feature_rows))):
        raise SubstrateError("arm A rna_feature_index is not 0..p-1 in row order")

    participants = [row["participant_id"] for row in participant_rows]
    if len(set(participants)) != len(participants):
        raise SubstrateError("arm A participant ids are not unique")
    folds = {row["participant_id"]: row["outer_fold"] for row in fold_rows}
    endpoints = {row["participant_id"]: row for row in endpoint_rows}
    if set(folds) != set(participants):
        raise SubstrateError("arm A fold roster differs from the molecular roster")
    if set(endpoints) != set(participants):
        raise SubstrateError("arm A endpoint roster differs from the molecular roster")
    disagreeing = [
        p for p in participants if endpoints[p]["outer_fold"] != folds[p]
    ]
    if disagreeing:
        raise SubstrateError(
            f"arm A outer_fold disagrees between the folds file and the endpoints "
            f"file for {len(disagreeing)} participants, first {disagreeing[:3]}"
        )

    for participant, record in endpoints.items():
        components = [int(record[name]) for name in NAS_COMPONENTS]
        if sum(components) != int(record["nash_crn_component_sum"]):
            raise SubstrateError(
                f"arm A {participant}: the three NASH-CRN components do not sum to "
                f"the deposited nash_crn_component_sum"
            )
        # Fail closed if lobular necrosis were ever folded into the NAS.
        if sum(components) + int(record[NOT_A_NAS_COMPONENT]) == int(
            record["nash_crn_component_sum"]
        ) and int(record[NOT_A_NAS_COMPONENT]) != 0:
            raise SubstrateError(
                f"arm A {participant}: lobular necrosis appears inside the NAS sum"
            )

    restricted = np.ascontiguousarray(values[:, shared_columns])
    np.save(output / "rna_shared_axis.npy", restricted)

    write_tsv(
        output / "participant_axis.tsv",
        participant_header + ["outer_fold"],
        [
            [row[name] for name in participant_header] + [folds[row["participant_id"]]]
            for row in participant_rows
        ],
    )
    endpoint_columns = [
        "participant_id",
        "outer_fold",
        *NAS_COMPONENTS,
        "nas_sum_deposited",
        NOT_A_NAS_COMPONENT,
        FOURTH_AXIS,
        "recorded_sex",
        "stage3",
        "source_stage5",
    ]
    write_tsv(
        output / "aspect_endpoints.tsv",
        endpoint_columns,
        [
            [
                participant,
                folds[participant],
                *[endpoints[participant][name] for name in NAS_COMPONENTS],
                endpoints[participant]["nash_crn_component_sum"],
                endpoints[participant][NOT_A_NAS_COMPONENT],
                endpoints[participant][FOURTH_AXIS],
                endpoints[participant]["recorded_sex"],
                endpoints[participant]["stage3"],
                endpoints[participant]["stage5"],
            ]
            for participant in participants
        ],
    )

    def census(name: str) -> dict[str, int]:
        return dict(
            sorted(Counter(endpoints[p][name] for p in participants).items())
        )

    receipt = {
        "schema_version": "masld-bench-showcase-arm-a-v1",
        "arm": "A",
        "role": "aspect_training_substrate",
        "source_series": "GSE267145",
        "biological_unit": "participant",
        "rows": len(participant_rows),
        "distinct_participants": len(set(participants)),
        "arithmetic_unit": "participant",
        "rna_source_samples": len({r["rna_source_sample_accession"] for r in participant_rows}),
        "h3k27ac_source_samples": len(
            {r["h3k27ac_source_sample_accession"] for r in participant_rows}
        ),
        "gene_universe": {
            "deposited_arm_a_features": len(feature_rows),
            "shared_with_arm_b": len(shared_genes),
            "filter_applied": "identity_intersection_with_arm_b_only",
            "expression_or_variance_filter_applied": False,
            "justification": (
                "The only filter is the identity intersection with the arm B "
                "axis. It reads no expression value and no outcome, so it cannot "
                "leak across an outer fold. An abundance or variance filter would "
                "have to be fit on all 99 participants at once, which the "
                "GSE267145 activation contract forbids: feature selection is fit "
                "inside outer training. The unrestricted 42,163-feature axis is "
                "left addressable for arm-A-only fits."
            ),
        },
        "aspects": {
            "nas_components": list(NAS_COMPONENTS),
            "nas_sum_reproduced_from_components": True,
            "lobular_necrosis_is_a_nas_component": False,
            "lobular_necrosis_carried_as": "separate deposited field",
            "fourth_axis": FOURTH_AXIS,
        },
        "census": {
            "steatosis": census("steatosis"),
            "ballooning": census("ballooning"),
            "lobular_inflammation": census("lobular_inflammation"),
            "lobular_necrosis": census(NOT_A_NAS_COMPONENT),
            "fibrosis": census(FOURTH_AXIS),
            "nas_sum": census("nash_crn_component_sum"),
            "recorded_sex": census("recorded_sex"),
            "outer_fold": dict(sorted(Counter(folds.values()).items())),
        },
        "matrix": {
            "shared_axis_shape": list(restricted.shape),
            "dtype": str(restricted.dtype),
            "normalization_applied": False,
        },
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    }
    write_json(output / "receipt.json", receipt)
    return receipt


def build_arm_b(
    *,
    source: Path,
    sra_metadata: Path,
    output: Path,
    shared_genes: Sequence[str],
    shared_columns: np.ndarray,
) -> dict[str, Any]:
    from scripts.freeze_gse135251_nash_source_inputs import read_rowname_offset_tsv

    molecular = source / "molecular"
    _, participant_rows = read_flat_tsv(molecular / "participant_axis.tsv")
    _, feature_rows = read_flat_tsv(molecular / "rna_feature_axis.tsv")
    _, outcome_rows = read_flat_tsv(source / "outcomes" / "participant_endpoints.tsv")
    _, fold_rows = read_flat_tsv(source / "folds" / "participant_outer_folds.tsv")
    values = np.load(molecular / "rna_values.npy")
    if values.shape != (len(participant_rows), len(feature_rows)):
        raise SubstrateError(
            f"arm B matrix {values.shape} does not match "
            f"{len(participant_rows)} participants x {len(feature_rows)} features"
        )
    participants = [row["participant_id"] for row in participant_rows]
    runs = [row["run_accession"] for row in participant_rows]
    if len(set(participants)) != len(participants) or len(set(runs)) != len(runs):
        raise SubstrateError("arm B participant or run accessions are not unique")

    outcomes = {row["participant_id"]: row for row in outcome_rows}
    folds = {row["participant_id"]: row for row in fold_rows}
    if set(outcomes) != set(participants) or set(folds) != set(participants):
        raise SubstrateError("arm B outcome or fold roster differs from the molecular roster")

    # Reconcile the 36 counts columns the frozen source carries no metadata for.
    # The bulk metadata table is R rowname-offset; the SRA table beside it is
    # flat.  Both invariants are asserted here rather than inherited.
    _, bulk_rows = read_rowname_offset_tsv(source / "raw" / "GSE135251_metadata.tsv")
    _, sra_rows = read_flat_tsv(sra_metadata)
    bulk_by_run = {row["sample_id"]: row for row in bulk_rows}
    sra_by_run = {row["sample_id"]: row for row in sra_rows}
    with (source / "raw" / "GSE135251_counts.tsv").open(encoding="utf-8") as handle:
        counts_header = handle.readline().rstrip("\n").split("\t")
    counts_runs = counts_header[1:]
    if len(set(counts_runs)) != len(counts_runs):
        raise SubstrateError("arm B counts columns are not unique")
    without_metadata = sorted(set(counts_runs) - set(bulk_by_run))
    if len(bulk_by_run) != len(participants) or set(bulk_by_run) != set(runs):
        raise SubstrateError("arm B bulk metadata roster differs from the frozen roster")
    unresolved = sorted(set(counts_runs) - set(sra_by_run))
    if unresolved:
        raise SubstrateError(
            f"{len(unresolved)} counts columns carry no record in either metadata "
            f"table, first {unresolved[:3]}; the 36 may not be dropped silently"
        )
    shared_fields = [f for f in sra_rows[0] if f in bulk_rows[0]]
    disagreements = {
        field: sum(
            1
            for run in bulk_by_run
            if bulk_by_run[run][field] != sra_by_run[run][field]
        )
        for field in shared_fields
    }
    disagreeing = {k: v for k, v in disagreements.items() if v}

    # Measure the contents of the column labelled SYMBOL.  Its name is not
    # evidence about what it holds.
    row_keys: list[str] = []
    with (source / "raw" / "GSE135251_counts.tsv").open(encoding="utf-8") as handle:
        handle.readline()
        for line in handle:
            if line.strip():
                row_keys.append(line.split("\t", 1)[0])
    symbol_column_contents = {
        "rows": len(row_keys),
        "unversioned_ensembl": sum(
            1 for k in row_keys if ENSEMBL_UNVERSIONED.fullmatch(k)
        ),
    }
    symbol_column_contents["not_ensembl"] = (
        symbol_column_contents["rows"] - symbol_column_contents["unversioned_ensembl"]
    )

    restricted = np.ascontiguousarray(values[:, shared_columns])
    np.save(output / "rna_shared_axis.npy", restricted)

    write_tsv(
        output / "participant_axis.tsv",
        ["participant_index", "participant_id", "run_accession",
         "outer_fold__nash_vs_nafl", "outer_fold__nash_vs_nafl_plus_control"],
        [
            [
                index,
                row["participant_id"],
                row["run_accession"],
                folds[row["participant_id"]]["outer_fold__nash_vs_nafl"],
                folds[row["participant_id"]]["outer_fold__nash_vs_nafl_plus_control"],
            ]
            for index, row in enumerate(participant_rows)
        ],
    )
    write_tsv(
        output / "observed_endpoints.tsv",
        ["participant_id", "run_accession", "nas_sum_observed", "fibrosis_stage",
         "group_in_paper", "disease"],
        [
            [
                row["participant_id"],
                row["run_accession"],
                outcomes[row["participant_id"]]["nas_score"],
                outcomes[row["participant_id"]]["fibrosis_stage"],
                outcomes[row["participant_id"]]["group_in_paper"],
                outcomes[row["participant_id"]]["disease"],
            ]
            for row in participant_rows
        ],
    )
    write_tsv(
        output / "counts_column_reconciliation.tsv",
        ["run_accession", "in_frozen_substrate", "metadata_provenance", "disease",
         "group_in_paper", "fibrosis_stage", "nas_score", "stage"],
        [
            [
                run,
                "yes" if run in bulk_by_run else "no",
                "bulk_and_sra" if run in bulk_by_run else "sra_only",
                sra_by_run[run]["disease"],
                sra_by_run[run]["group_in_paper"],
                sra_by_run[run]["Fibrosis_stage"],
                sra_by_run[run]["nas_score"],
                sra_by_run[run]["Stage"],
            ]
            for run in counts_runs
        ],
    )

    def sra_census(runs_subset: Sequence[str], field: str) -> dict[str, int]:
        return dict(sorted(Counter(sra_by_run[r][field] for r in runs_subset).items()))

    receipt = {
        "schema_version": "masld-bench-showcase-arm-b-v1",
        "arm": "B",
        "role": "external_check_substrate",
        "source_series": "GSE135251",
        "reused_frozen_source": True,
        "source_rebuilt": False,
        "biological_unit": "participant",
        "rows": len(participant_rows),
        "distinct_participants": len(set(participants)),
        "distinct_run_accessions": len(set(runs)),
        "arithmetic_unit": "participant",
        "aspects": {
            "components_deposited": False,
            "nas_sum_deposited": True,
            "fibrosis_stage_deposited": True,
            "note": (
                "GSE135251 deposits the NAS sum and a fibrosis stage but not the "
                "three components. The external check is therefore a "
                "reconstruction check: three predicted aspects summed against an "
                "observed sum, plus the fourth axis scored directly."
            ),
        },
        "counts_column_reconciliation": {
            "counts_columns_deposited": len(counts_runs),
            "columns_in_frozen_substrate": len(participants),
            "columns_without_bulk_metadata": len(without_metadata),
            "columns_without_any_metadata": len(unresolved),
            "runs_without_bulk_metadata": without_metadata,
            "resolved_against": str(sra_metadata),
            "resolved_against_sha256": sha256_file(sra_metadata),
            "bulk_metadata_is_a_strict_subset_of_the_sra_table": set(bulk_by_run)
            <= set(sra_by_run),
            "shared_fields_compared": len(shared_fields),
            "shared_field_disagreements_on_the_180": disagreeing,
            "finding": (
                "The 36 are not unlabelled at the source. Every one of them "
                "carries disease, group_in_paper, Fibrosis_stage, nas_score and "
                "Stage in the 216-row SRA table, and that table agrees with the "
                "bulk table on all 32 fields of all 180 shared runs. The frozen "
                "n=180 is a subset of an available n=216 because the freeze read "
                "the 180-row bulk metadata file, not because the source is "
                "missing labels."
            ),
            "action_taken": (
                "None. The frozen source is reused unchanged as instructed. The "
                "36 are recorded here with their labels as a documented, "
                "recoverable extension and are never silently pooled into the "
                "primary roster."
            ),
            "census_of_the_36": {
                field: sra_census(without_metadata, field)
                for field in ("disease", "group_in_paper", "Fibrosis_stage",
                              "nas_score", "Stage")
            },
            "census_of_all_216": {
                field: sra_census(counts_runs, field)
                for field in ("disease", "group_in_paper", "Fibrosis_stage",
                              "nas_score")
            },
        },
        "symbol_column_contents": symbol_column_contents,
        "gene_universe": {
            "deposited_arm_b_features": len(feature_rows),
            "shared_with_arm_a": len(shared_genes),
        },
        "census": {
            "nas_sum_observed": dict(
                sorted(Counter(outcomes[p]["nas_score"] for p in participants).items())
            ),
            "fibrosis_stage": dict(
                sorted(
                    Counter(outcomes[p]["fibrosis_stage"] for p in participants).items()
                )
            ),
            "group_in_paper": dict(
                sorted(
                    Counter(outcomes[p]["group_in_paper"] for p in participants).items()
                )
            ),
        },
        "matrix": {
            "shared_axis_shape": list(restricted.shape),
            "dtype": str(restricted.dtype),
            "normalization_applied": False,
        },
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    }
    write_json(output / "receipt.json", receipt)
    return receipt


def build_arm_c(
    *,
    series_matrix: Path,
    run_info: Path,
    single_nucleus: Path,
    donor_pairing: Path,
    output: Path,
) -> dict[str, Any]:
    accessions, fields = parse_series_matrix_characteristics(series_matrix)
    missing_key = [a for a in accessions if "patient_id" not in fields[a]]
    if missing_key:
        raise SubstrateError(
            f"{len(missing_key)} GSMs carry no patient id, first {missing_key[:3]}"
        )
    absent_saf = [a for a in accessions if "saf_score" not in fields[a]]
    if absent_saf:
        raise SubstrateError(
            f"the saf score field is absent for {len(absent_saf)} GSMs; the frozen "
            f"spec records it present on every sample, so this parse is wrong"
        )

    graded: dict[str, tuple[int, int, int]] = {}
    ungraded: dict[str, str] = {}
    for accession in accessions:
        match = SAF_TRIPLE.fullmatch(fields[accession]["saf_score"])
        if match:
            graded[accession] = tuple(int(v) for v in match.groups())  # type: ignore[assignment]
        else:
            ungraded[accession] = fields[accession]["saf_score"]

    donor_of_gsm = {a: fields[a]["patient_id"] for a in accessions}
    donors = sorted(set(donor_of_gsm.values()))
    donor_status: dict[str, set[str]] = defaultdict(set)
    for accession in accessions:
        donor_status[donor_of_gsm[accession]].add(
            fields[accession].get("disease_status", "")
        )
    split_status = {d: sorted(v) for d, v in donor_status.items() if len(v) > 1}
    if split_status:
        raise SubstrateError(f"donors carry conflicting disease status: {split_status}")
    donor_triples: dict[str, set[tuple[int, int, int]]] = defaultdict(set)
    for accession, triple in graded.items():
        donor_triples[donor_of_gsm[accession]].add(triple)
    conflicted = {d: sorted(v) for d, v in donor_triples.items() if len(v) > 1}
    if conflicted:
        raise SubstrateError(f"within-donor SAF conflicts: {conflicted}")
    graded_donors = sorted(donor_triples)
    ungraded_donors = sorted(set(donors) - set(graded_donors))

    runs = read_sra_run_info(run_info)
    run_to_gsm = {row["run_accession"]: row["sample_accession"] for row in runs}
    unknown = sorted(set(run_to_gsm.values()) - set(accessions))
    if unknown:
        raise SubstrateError(f"run table names GSMs absent from the series: {unknown[:3]}")

    obs = read_h5ad_obs_columns(single_nucleus, ["cell_type", "sample"])
    nucleus_runs = obs["sample"]
    nucleus_lineages = obs["cell_type"]
    orphan = sorted(set(nucleus_runs) - set(run_to_gsm))
    if orphan:
        raise SubstrateError(
            f"{len(orphan)} single-nucleus run accessions are absent from the "
            f"primary run table, first {orphan[:3]}"
        )
    lineage_axis = sorted(set(nucleus_lineages))

    counts: dict[str, Counter[str]] = defaultdict(Counter)
    donor_runs: dict[str, set[str]] = defaultdict(set)
    for run, lineage in zip(nucleus_runs, nucleus_lineages):
        donor = donor_of_gsm[run_to_gsm[run]]
        counts[donor][lineage] += 1
        donor_runs[donor].add(run)
    donor_gsms: dict[str, set[str]] = defaultdict(set)
    for accession, donor in donor_of_gsm.items():
        donor_gsms[donor].add(accession)

    # Cross-check the derived map against the project's donor pairing table.
    # The pairing table is a derived output file, so it is compared, never adopted.
    _, pairing_rows = read_flat_tsv_or_csv(donor_pairing)
    paired: dict[str, str] = {}
    for row in pairing_rows:
        for run in str(row["rna_srrs"]).split(";"):
            run = run.strip()
            if run:
                paired[run] = row["donor_id"]
    derived = {run: donor_of_gsm[gsm] for run, gsm in run_to_gsm.items()}
    pairing_check = {
        "pairing_table_rows": len(pairing_rows),
        "pairing_table_distinct_donors": len({r["donor_id"] for r in pairing_rows}),
        "pairing_table_runs": len(paired),
        "derived_runs": len(derived),
        "runs_only_in_derived": sorted(set(derived) - set(paired)),
        "runs_only_in_pairing_table": sorted(set(paired) - set(derived)),
        "donor_disagreements": {
            run: {"derived": derived[run], "pairing_table": paired[run]}
            for run in sorted(set(derived) & set(paired))
            if derived[run] != paired[run]
        },
        "adopted": False,
        "role": "cross_check_only",
    }

    write_tsv(
        output / "run_to_sample.tsv",
        ["run_accession", "experiment_accession", "sample_accession", "donor_id",
         "in_single_nucleus_matrix"],
        [
            [
                row["run_accession"],
                row["experiment_accession"],
                row["sample_accession"],
                donor_of_gsm[row["sample_accession"]],
                "yes" if row["run_accession"] in set(nucleus_runs) else "no",
            ]
            for row in sorted(runs, key=lambda r: r["run_accession"])
        ],
    )
    write_tsv(
        output / "sample_to_donor.tsv",
        ["sample_accession", "donor_id", "disease_status", "saf_score_literal",
         "saf_grade_state"],
        [
            [
                accession,
                donor_of_gsm[accession],
                fields[accession].get("disease_status", ""),
                fields[accession]["saf_score"],
                "saf_graded" if accession in graded else "ungraded",
            ]
            for accession in accessions
        ],
    )
    write_tsv(
        output / "donor_aspects.tsv",
        ["donor_id", "saf_grade_state", "disease_status", "saf_steatosis_S",
         "saf_activity_A", "saf_fibrosis_F", "nas_equivalent_S_plus_A",
         "n_samples", "n_runs", "n_runs_in_matrix", "n_nuclei"],
        [
            [
                donor,
                "saf_graded" if donor in donor_triples else "ungraded",
                next(iter(donor_status[donor])),
                *(
                    list(next(iter(donor_triples[donor])))
                    if donor in donor_triples
                    else ["", "", ""]
                ),
                (
                    sum(next(iter(donor_triples[donor]))[:2])
                    if donor in donor_triples
                    else ""
                ),
                len(donor_gsms[donor]),
                len([r for r, g in run_to_gsm.items() if g in donor_gsms[donor]]),
                len(donor_runs.get(donor, ())),
                sum(counts.get(donor, Counter()).values()),
            ]
            for donor in donors
        ],
    )
    write_tsv(
        output / "donor_lineage_counts.tsv",
        ["donor_id", "n_nuclei", *lineage_axis],
        [
            [donor, sum(counts[donor].values()),
             *[counts[donor][lineage] for lineage in lineage_axis]]
            for donor in sorted(counts)
        ],
    )
    write_tsv(
        output / "donor_lineage_fractions.tsv",
        ["donor_id", "n_nuclei", *lineage_axis],
        [
            [
                donor,
                sum(counts[donor].values()),
                *[
                    f"{counts[donor][lineage] / sum(counts[donor].values()):.10f}"
                    for lineage in lineage_axis
                ],
            ]
            for donor in sorted(counts)
        ],
    )

    graded_with_nuclei = sorted(set(graded_donors) & set(counts))
    triples = [next(iter(donor_triples[d])) for d in graded_donors]
    receipt = {
        "schema_version": "masld-bench-showcase-arm-c-v1",
        "arm": "C",
        "role": "aspect_and_lineage_cross_check",
        "is_training_data": False,
        "source_series": "GSE202379",
        "biological_unit": "donor",
        "unit_census": {
            "sample_rows_gsm": len(accessions),
            "sequencing_runs": len(runs),
            "runs_in_single_nucleus_matrix": len(set(nucleus_runs)),
            "runs_absent_from_the_matrix": sorted(
                set(run_to_gsm) - set(nucleus_runs)
            ),
            "distinct_donors": len(donors),
            "saf_graded_donors": len(graded_donors),
            "ungraded_donors": len(ungraded_donors),
            "graded_donors_with_nuclei": len(graded_with_nuclei),
            "nuclei": int(len(nucleus_runs)),
            "arithmetic_unit": "donor",
            "statement": (
                f"{len(donors)} donors across {len(accessions)} GSMs and "
                f"{len(runs)} runs; {len(graded_donors)} are SAF-graded. Sample "
                f"rows are never the inferential unit."
            ),
        },
        "stage_source": {
            "used": "primary GEO series matrix saf score field",
            "series_matrix_sha256": sha256_file(series_matrix),
            "donor_fstage_documented_tsv_used": False,
            "why_barred": (
                "It invents seven stages, overrides one and drops one against the "
                "primary SAF, and is keyed on runs rather than donors."
            ),
            "ungraded_definition": "fails ^S\\d+A\\d+F\\d+$, never an absent field",
            "ungraded_literals": dict(sorted(Counter(ungraded.values()).items())),
            "ungraded_donors": ungraded_donors,
        },
        "aspect_resolution": {
            "steatosis_S": "0-3, directly observed",
            "activity_A": "0-4, directly observed",
            "fibrosis_F": "0-4, directly observed",
            "ballooning_and_lobular_inflammation_separable": False,
            "limitation": (
                "SAF deposits activity A as one number. A is the ballooning plus "
                "lobular inflammation sum, so this arm can cross-check steatosis "
                "and the ballooning-plus-inflammation SUM, but cannot resolve "
                "ballooning against lobular inflammation individually."
            ),
            "nas_equivalence": (
                "S + A occupies the same slot as the NASH-CRN NAS because SAF's A "
                "is defined as ballooning (0-2) plus lobular inflammation (0-2). "
                "SAF and NASH-CRN remain different published scoring systems and "
                "this correspondence is a scale alignment, not a claim that the "
                "two were read by the same protocol."
            ),
        },
        "lineage": {
            "axis": lineage_axis,
            "levels": len(lineage_axis),
            "directly_observed": True,
            "operator": "single-nucleus cell_type annotation",
            "is_a_deconvolution_estimate": False,
            "note": (
                "This is the only arm where aspect and lineage are both directly "
                "observed. Its lineage axis has "
                f"{len(lineage_axis)} levels and is NOT the 16-level REF_HUMAN "
                "BayesPrism axis; the two are not interchangeable."
            ),
        },
        "donor_saf_census": {
            "S": dict(sorted(Counter(t[0] for t in triples).items())),
            "A": dict(sorted(Counter(t[1] for t in triples).items())),
            "F": dict(sorted(Counter(t[2] for t in triples).items())),
            "S_plus_A": dict(sorted(Counter(t[0] + t[1] for t in triples).items())),
        },
        "donor_pairing_cross_check": pairing_check,
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    }
    write_json(output / "receipt.json", receipt)
    return receipt


def read_flat_tsv_or_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter=delimiter))
    if len(rows) < 2:
        raise SubstrateError(f"table has no data rows: {path}")
    header = rows[0]
    for index, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise SubstrateError(
                f"{path} line {index} carries {len(row)} fields, not {len(header)}"
            )
    return header, [dict(zip(header, row)) for row in rows[1:]]


def read_sha256_manifest(path: Path) -> dict[Path, str]:
    manifest: dict[Path, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, _, name = line.partition("  ")
        manifest[Path(name)] = digest
    if not manifest:
        raise SubstrateError(f"empty manifest: {path}")
    return manifest


def build_lineage_substrate(
    *, results_root: Path, manifest_path: Path, output: Path
) -> dict[str, Any]:
    tables = sorted(results_root.glob("*/*_bayesprism_proportions.tsv"))
    if not tables:
        raise SubstrateError(f"no BayesPrism proportions under {results_root}")
    # The glob is pinned to a manifest so a table appearing or vanishing between
    # the pin and the run is an error rather than a silently different substrate.
    manifest = read_sha256_manifest(manifest_path)
    resolved = {table.resolve(): table for table in tables}
    unpinned = sorted(str(p) for p in set(resolved) - set(manifest))
    missing = sorted(str(p) for p in set(manifest) - set(resolved))
    if unpinned or missing:
        raise SubstrateError(
            f"BayesPrism roster drifted from the pinned manifest: "
            f"unpinned={unpinned}, missing={missing}"
        )
    admitted: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for table in tables:
        lineages, keys, matrix = read_bayesprism_proportions(table)
        sums = matrix.sum(axis=1)
        record = {
            "cohort": table.parent.name,
            "file": table.name,
            "source_path": str(table),
            "source_sha256": sha256_file(table),
            "pinned_sha256_matches": sha256_file(table)
            == manifest[table.resolve()],
            "rows": len(keys),
            "distinct_keys": len(set(keys)),
            "lineage_levels": len(lineages),
            "row_sum_min": float(sums.min()),
            "row_sum_max": float(sums.max()),
            "row_sums_are_one": bool(np.allclose(sums, 1.0, atol=1e-9, rtol=0.0)),
            "exact_zero_cells": int((matrix == 0.0).sum()),
            "key_looks_like_a_run_accession": bool(
                all(RUN_ACCESSION.fullmatch(k) for k in keys)
            ),
            "key_looks_like_a_sample_accession": bool(
                all(SAMPLE_ACCESSION.fullmatch(k) for k in keys)
            ),
        }
        if not record["pinned_sha256_matches"]:
            raise SubstrateError(f"{table}: content differs from the pinned manifest")
        if not record["row_sums_are_one"]:
            raise SubstrateError(
                f"{table}: row sums are not 1 "
                f"({record['row_sum_min']}, {record['row_sum_max']})"
            )
        if tuple(lineages) == REF_HUMAN_LINEAGES:
            destination = output / "ref_human_16_lineage" / table.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(table, destination)
            record["frozen_copy"] = f"ref_human_16_lineage/{table.name}"
            record["frozen_copy_sha256"] = sha256_file(destination)
            if record["frozen_copy_sha256"] != record["source_sha256"]:
                raise SubstrateError(f"{table}: frozen copy differs from the source")
            admitted.append(record)
        else:
            record["excluded_because"] = (
                f"lineage axis has {len(lineages)} levels and is not the REF_HUMAN "
                f"16-level axis"
            )
            record["lineage_axis"] = lineages
            excluded.append(record)
    write_tsv(
        output / "lineage_axis.tsv",
        ["lineage_index", "lineage"],
        [[i, name] for i, name in enumerate(REF_HUMAN_LINEAGES)],
    )
    receipt = {
        "schema_version": "masld-bench-showcase-lineage-v1",
        "role": "lineage_substrate",
        "operator": "BayesPrism",
        "music_used": False,
        "why_not_music": (
            "MuSiC is recorded degenerate on real bulk in this project: "
            "cholangiocytes exactly zero in 97.8% of samples under both "
            "references."
        ),
        "ref_human_axis": list(REF_HUMAN_LINEAGES),
        "ref_human_levels": len(REF_HUMAN_LINEAGES),
        "cohorts_frozen": len(admitted),
        "cohorts_excluded": len(excluded),
        "admitted": admitted,
        "excluded": excluded,
        "row_sums_verified_to_one": True,
        "row_sum_tolerance": 1e-9,
        "validated_for": "composition_only",
        "not_validated_for": "between_lineage_statements",
        "subcomposition_caveat": (
            "The operator is validated for composition only. The W2 macrophage "
            "result reversed sign, +0.225 to -0.165, under subcomposition, so any "
            "between-lineage statement needs its own subcomposition check before "
            "it is made."
        ),
        "shared_reference_caveat": (
            "Every admitted cohort was deconvolved against one shared REF_HUMAN "
            "reference. Cross-cohort agreement therefore shows one operator "
            "transferring, not independent measurements agreeing."
        ),
        "donor_key_caveat": (
            "Every deconvolution cohort here is keyed by sequencing run or GEO "
            "sample, with no donor key. No donor count may be claimed from these "
            "tables. rectangle_comparison/pseudobulk/donor_folds.tsv already "
            "mislabels SRR accessions as a donor column."
        ),
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    }
    write_json(output / "receipt.json", receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm-a-fixture", type=Path, required=True)
    parser.add_argument("--arm-a-activation", type=Path, required=True)
    parser.add_argument("--arm-b-source", type=Path, required=True)
    parser.add_argument("--arm-b-sra-metadata", type=Path, required=True)
    parser.add_argument("--arm-c-series-matrix", type=Path, required=True)
    parser.add_argument("--arm-c-run-info", type=Path, required=True)
    parser.add_argument("--arm-c-single-nucleus", type=Path, required=True)
    parser.add_argument("--arm-c-donor-pairing", type=Path, required=True)
    parser.add_argument("--arm-c-frozen-spec", type=Path, required=True)
    parser.add_argument("--bayesprism-results", type=Path, required=True)
    parser.add_argument("--bayesprism-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = args.output
    for name in ("arm_a", "arm_b", "arm_c", "lineage", "shared", "sources"):
        (output / name).mkdir(parents=True, exist_ok=True)

    _, arm_a_features = read_flat_tsv(
        args.arm_a_fixture / "molecular" / "rna_feature_axis.tsv"
    )
    _, arm_b_features = read_flat_tsv(
        args.arm_b_source / "molecular" / "rna_feature_axis.tsv"
    )
    shared, a_columns, b_columns = build_shared_gene_axis(
        [row["stable_gene_id"] for row in arm_a_features],
        [row["stable_gene_id"] for row in arm_b_features],
    )
    write_tsv(
        output / "shared" / "shared_gene_axis.tsv",
        ["shared_index", "stable_gene_id", "arm_a_source_index", "arm_b_source_index"],
        [
            [index, gene, int(a_columns[index]), int(b_columns[index])]
            for index, gene in enumerate(shared)
        ],
    )

    arm_a = build_arm_a(
        fixture=args.arm_a_fixture,
        activation=args.arm_a_activation,
        output=output / "arm_a",
        shared_genes=shared,
        shared_columns=a_columns,
    )
    arm_b = build_arm_b(
        source=args.arm_b_source,
        sra_metadata=args.arm_b_sra_metadata,
        output=output / "arm_b",
        shared_genes=shared,
        shared_columns=b_columns,
    )
    arm_c = build_arm_c(
        series_matrix=args.arm_c_series_matrix,
        run_info=args.arm_c_run_info,
        single_nucleus=args.arm_c_single_nucleus,
        donor_pairing=args.arm_c_donor_pairing,
        output=output / "arm_c",
    )
    lineage = build_lineage_substrate(
        results_root=args.bayesprism_results,
        manifest_path=args.bayesprism_manifest,
        output=output / "lineage",
    )

    # Arm C must reproduce the fibrosis grades of the pre-registered frozen spec
    # independently.  A parse that disagrees is a parse that is wrong.
    _, spec_rows = read_flat_tsv(args.arm_c_frozen_spec / "spec" / "graded_donors.tsv")
    spec_stage = {
        row["donor_id"]: row["saf_fibrosis_stage"]
        for row in spec_rows
        if row["f_stage_source"] == "saf_graded"
    }
    _, derived_rows = read_flat_tsv(output / "arm_c" / "donor_aspects.tsv")
    derived_stage = {
        row["donor_id"]: row["saf_fibrosis_F"]
        for row in derived_rows
        if row["saf_grade_state"] == "saf_graded"
    }
    if spec_stage != derived_stage:
        raise SubstrateError(
            "arm C fibrosis grades do not reproduce the pre-registered frozen spec"
        )

    # Freeze hashed read-only copies of every external input that is small
    # enough to carry.  The single-nucleus matrix is recorded by digest only.
    external = {
        "gse202379_series_matrix.txt.gz": args.arm_c_series_matrix,
        "gse202379_SraRunInfo.csv": args.arm_c_run_info,
        "gse202379_donor_pairing.csv": args.arm_c_donor_pairing,
        "gse135251_sra_metadata.tsv": args.arm_b_sra_metadata,
        "bayesprism_sources.sha256": args.bayesprism_manifest,
    }
    frozen_sources = {}
    for name, path in external.items():
        shutil.copyfile(path, output / "sources" / name)
        digest = sha256_file(output / "sources" / name)
        if digest != sha256_file(path):
            raise SubstrateError(f"frozen copy of {path} differs from the source")
        frozen_sources[name] = {"source_path": str(path), "sha256": digest}
    frozen_sources["gse202379_single_nucleus.h5ad"] = {
        "source_path": str(args.arm_c_single_nucleus),
        "sha256": sha256_file(args.arm_c_single_nucleus),
        "copied": False,
        "why_not_copied": "693 MB; recorded by digest and read for obs columns only",
    }
    write_json(output / "sources" / "external_inputs.json", frozen_sources)

    summary = {
        "schema_version": "masld-bench-showcase-substrate-v1",
        "purpose": "aspect_and_lineage_resolved_gene_model_substrate",
        "prepared_only": True,
        "model_fitted": False,
        "metrics_calculated": False,
        "shared_gene_axis": {
            "size": len(shared),
            "derivation": (
                "Inner join on unversioned Ensembl stable gene id between the "
                f"deposited arm A axis ({len(arm_a_features)}) and the frozen arm "
                f"B bridge ({len(arm_b_features)}). Both axes are already "
                "unversioned and duplicate-free, so the join needs no collapse "
                "rule. Identity only: no expression value and no outcome is read."
            ),
            "arm_a_axis": len(arm_a_features),
            "arm_b_axis": len(arm_b_features),
            "arm_a_only": len(arm_a_features) - len(shared),
            "arm_b_only": len(arm_b_features) - len(shared),
        },
        "realised_units": {
            "arm_a_participants": arm_a["distinct_participants"],
            "arm_b_participants": arm_b["distinct_participants"],
            "arm_c_donors": arm_c["unit_census"]["distinct_donors"],
            "arm_c_saf_graded_donors": arm_c["unit_census"]["saf_graded_donors"],
            "lineage_cohorts_frozen": lineage["cohorts_frozen"],
        },
        "status": "passed",
    }
    write_json(output / "substrate_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
