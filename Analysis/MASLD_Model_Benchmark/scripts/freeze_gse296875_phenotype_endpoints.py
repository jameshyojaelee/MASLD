#!/usr/bin/env python3
"""Freeze GSE296875 donor endpoint rows with explicit typed missingness masks.

The deposited sheet writes an unobserved donor field as the literal string
``NA``.  Every field therefore has a full column of non-blank cells, and a
blank-based mask would silently adopt missing pathology to an observed
score.  Each endpoint here carries its own ``observed`` flag; a donor who is
not observed carries a null value and never a zero, a stage, or a healthy
default.

Values stay on their source-native scale.  Fibrosis keeps the deposited
categorical levels and is mapped to a binary any-fibrosis endpoint only through
the exact level lists frozen in the TaskSpec; an unlisted level is a hard
failure rather than a silent negative.  No ordinal fibrosis stage is inferred.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import tomllib
from typing import Any

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


MISSING_SENTINELS = ("NA",)
DONOR_FIELD = "sample_id"
COVARIATE_FIELDS = (
    "age_in_yr",
    "reported_sex",
    "reported_race",
    "BMI",
    "height_m",
    "weight_kg",
    "cause_of_death",
    "genetic_similarity",
)
NUMERIC_FIELDS = frozenset(
    {"age_in_yr", "BMI", "height_m", "weight_kg", "steatosis_numeric", "n_cells"}
)


class PhenotypeEndpointError(ValueError):
    """Raised when the supplement contradicts the frozen phenotype requirements."""


def _clean(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def read_sheet(cells_path: Path) -> tuple[list[str], list[dict[str, str | None]]]:
    """Return the header and one dict per data row of the frozen sheet dump."""

    rows = json.loads(cells_path.read_text(encoding="utf-8"))
    header_index = next(
        (
            index
            for index, row in enumerate(rows)
            if sum(1 for cell in row if _clean(cell["value"]) is not None) >= 2
        ),
        None,
    )
    if header_index is None:
        raise PhenotypeEndpointError("frozen sheet dump has no header row")
    header = [_clean(cell["value"]) or "" for cell in rows[header_index]]
    while header and header[-1] == "":
        header.pop()
    if len(set(header)) != len(header):
        raise PhenotypeEndpointError("frozen sheet header is not unique")
    records = []
    for row in rows[header_index + 1 :]:
        values = [_clean(cell["value"]) for cell in row]
        if not any(value is not None for value in values):
            continue
        values.extend([None] * (len(header) - len(values)))
        records.append(dict(zip(header, values[: len(header)])))
    return header, records


def _observed(value: str | None) -> bool:
    return value is not None and value not in MISSING_SENTINELS


def _numeric(field: str, value: str) -> float:
    try:
        return float(value)
    except ValueError as error:
        raise PhenotypeEndpointError(f"{field} is not numeric: {value!r}") from error


def build_endpoints(
    records: list[dict[str, str | None]],
    roster: list[str],
    spec: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Attach typed endpoints and masks to exactly the analyzed donor roster."""

    by_donor: dict[str, dict[str, str | None]] = {}
    for record in records:
        donor = record.get(DONOR_FIELD)
        if donor is None:
            raise PhenotypeEndpointError("a supplement row carries no donor identifier")
        if donor in by_donor:
            raise PhenotypeEndpointError(f"donor appears twice in the supplement: {donor}")
        by_donor[donor] = record
    absent = [donor for donor in roster if donor not in by_donor]
    if absent:
        raise PhenotypeEndpointError(f"analyzed donors absent from the supplement: {absent}")
    deposited_not_analyzed = sorted(set(by_donor).difference(roster))

    negatives = set(spec["fibrosis"]["negative_source_values"])
    positives = set(spec["fibrosis"]["positive_source_values"])
    if negatives & positives:
        raise PhenotypeEndpointError("fibrosis level lists overlap")

    endpoints = []
    for donor in roster:
        record = by_donor[donor]
        row: dict[str, Any] = {"donor_id": donor, "well_id": record.get("well")}

        steatosis = record.get("steatosis_numeric")
        steatosis_observed = _observed(steatosis)
        row["steatosis_observed"] = steatosis_observed
        row["steatosis_numeric"] = (
            _numeric("steatosis_numeric", steatosis) if steatosis_observed else None
        )
        row["steatosis_categorical"] = (
            record.get("steatosis_categorical")
            if _observed(record.get("steatosis_categorical"))
            else None
        )
        row["steatosis_status"] = (
            record.get("steatosis_status")
            if _observed(record.get("steatosis_status"))
            else None
        )

        fibrosis = record.get("fibrosis_categorical")
        fibrosis_observed = _observed(fibrosis)
        row["fibrosis_observed"] = fibrosis_observed
        row["fibrosis_categorical"] = fibrosis if fibrosis_observed else None
        if fibrosis_observed:
            if fibrosis in positives:
                row["fibrosis_any"] = True
            elif fibrosis in negatives:
                row["fibrosis_any"] = False
            else:
                raise PhenotypeEndpointError(
                    f"fibrosis level is in neither frozen list: {fibrosis!r}"
                )
        else:
            row["fibrosis_any"] = None
        row["fibrosis_status"] = (
            record.get("fibrosis_status")
            if _observed(record.get("fibrosis_status"))
            else None
        )

        for field in COVARIATE_FIELDS:
            raw = record.get(field)
            observed = _observed(raw)
            row[f"{field}_observed"] = observed
            if not observed:
                row[field] = None
            elif field in NUMERIC_FIELDS:
                row[field] = _numeric(field, raw)
            else:
                row[field] = raw
        if row["age_in_yr_observed"]:
            row["is_adult"] = row["age_in_yr"] >= 18
        else:
            row["is_adult"] = None
        endpoints.append(row)

    steatosis_observed = [row for row in endpoints if row["steatosis_observed"]]
    fibrosis_observed = [row for row in endpoints if row["fibrosis_observed"]]
    ages = [row["age_in_yr"] for row in endpoints if row["age_in_yr_observed"]]
    summary = {
        "analyzed_donors": len(endpoints),
        "deposited_donors_not_analyzed": deposited_not_analyzed,
        "steatosis_observed_donors": len(steatosis_observed),
        "steatosis_missing_donors": len(endpoints) - len(steatosis_observed),
        "steatosis_missing_donor_ids": [
            row["donor_id"] for row in endpoints if not row["steatosis_observed"]
        ],
        "fibrosis_observed_donors": len(fibrosis_observed),
        "fibrosis_missing_donors": len(endpoints) - len(fibrosis_observed),
        "fibrosis_missing_donor_ids": [
            row["donor_id"] for row in endpoints if not row["fibrosis_observed"]
        ],
        "fibrosis_positive_donors": sum(1 for row in fibrosis_observed if row["fibrosis_any"]),
        "fibrosis_negative_donors": sum(
            1 for row in fibrosis_observed if row["fibrosis_any"] is False
        ),
        "fibrosis_observed_levels": sorted(
            {row["fibrosis_categorical"] for row in fibrosis_observed}
        ),
        "steatosis_observed_levels": sorted(
            {
                row["steatosis_categorical"]
                for row in endpoints
                if row["steatosis_categorical"] is not None
            }
        ),
        "steatosis_numeric_range": (
            [
                min(row["steatosis_numeric"] for row in steatosis_observed),
                max(row["steatosis_numeric"] for row in steatosis_observed),
            ]
            if steatosis_observed
            else None
        ),
        "covariate_observed_donors": {
            field: sum(1 for row in endpoints if row[f"{field}_observed"])
            for field in COVARIATE_FIELDS
        },
        "reported_sex_counts": _counts(endpoints, "reported_sex"),
        "reported_race_counts": _counts(endpoints, "reported_race"),
        "genetic_similarity_counts": _counts(endpoints, "genetic_similarity"),
        "minimum_age_years": min(ages) if ages else None,
        "maximum_age_years": max(ages) if ages else None,
        "donors_under_18": sum(1 for row in endpoints if row["is_adult"] is False),
        "donors_under_18_ids": [row["donor_id"] for row in endpoints if row["is_adult"] is False],
        "adult_donors": sum(1 for row in endpoints if row["is_adult"] is True),
        "missing_sentinels": list(MISSING_SENTINELS),
        "recoding": spec["steatosis"]["recoding"],
        "ordinal_stage_inference": spec["fibrosis"]["ordinal_stage_inference"],
    }
    return endpoints, summary


def _counts(endpoints: list[dict[str, Any]], field: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in endpoints:
        if row[f"{field}_observed"]:
            counts[row[field]] = counts.get(row[field], 0) + 1
    return dict(sorted(counts.items()))


def _write_tsv(path: Path, header: list[str], endpoints: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in endpoints:
            fields = []
            for name in header:
                value = row[name]
                if value is None:
                    fields.append("")
                elif isinstance(value, bool):
                    fields.append("true" if value else "false")
                elif isinstance(value, float) and value.is_integer():
                    fields.append(str(int(value)))
                else:
                    fields.append(str(value))
            writer.writerow(fields)
    path.chmod(0o440)


def run(*, cells: Path, join_lock: Path, taskspec: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise PhenotypeEndpointError("refusing to overwrite phenotype-endpoint artifact")
    with taskspec.open("rb") as handle:
        spec = tomllib.load(handle)
    with (join_lock / "donor_census.tsv").open(encoding="utf-8", newline="") as handle:
        census = list(csv.DictReader(handle, delimiter="\t"))
    roster = [row["donor_id"] for row in census]
    join_well = {row["donor_id"]: row["well_id"] for row in census}
    if len(set(roster)) != len(roster):
        raise PhenotypeEndpointError("the frozen donor census repeats a donor")
    if len(roster) != spec["analyzed_donors"]:
        raise PhenotypeEndpointError(
            f"donor roster differs: join {len(roster)}, TaskSpec {spec['analyzed_donors']}"
        )

    header, records = read_sheet(cells)
    for field in spec["source_fields"]:
        if field not in header:
            raise PhenotypeEndpointError(f"declared source field is absent: {field}")
    endpoints, summary = build_endpoints(records, roster, spec)

    disagreeing = [
        row["donor_id"] for row in endpoints if row["well_id"] != join_well[row["donor_id"]]
    ]
    if disagreeing:
        raise PhenotypeEndpointError(
            f"supplement well disagrees with the frozen join for donors: {disagreeing}"
        )
    summary["supplement_well_agrees_with_frozen_join"] = True

    for endpoint, observed_key, missing_key in (
        ("steatosis", "steatosis_observed_donors", "steatosis_missing_donors"),
        ("fibrosis", "fibrosis_observed_donors", "fibrosis_missing_donors"),
    ):
        if summary[observed_key] != spec[endpoint]["observed_donors"]:
            raise PhenotypeEndpointError(
                f"{endpoint} observed donors differ: observed {summary[observed_key]}, "
                f"frozen {spec[endpoint]['observed_donors']}"
            )
        if summary[missing_key] != spec[endpoint]["missing_donors"]:
            raise PhenotypeEndpointError(
                f"{endpoint} missing donors differ: observed {summary[missing_key]}, "
                f"frozen {spec[endpoint]['missing_donors']}"
            )
    if summary["donors_under_18"] != spec["donors_under_18"]:
        raise PhenotypeEndpointError(
            f"donors under 18 differ: observed {summary['donors_under_18']}, "
            f"frozen {spec['donors_under_18']}"
        )
    if summary["minimum_age_years"] != spec["minimum_age_years"]:
        raise PhenotypeEndpointError("minimum age differs from the frozen TaskSpec")
    if summary["maximum_age_years"] != spec["maximum_age_years"]:
        raise PhenotypeEndpointError("maximum age differs from the frozen TaskSpec")
    unlisted = set(summary["fibrosis_observed_levels"]).difference(
        set(spec["fibrosis"]["negative_source_values"])
        | set(spec["fibrosis"]["positive_source_values"])
    )
    if unlisted:
        raise PhenotypeEndpointError(f"fibrosis levels outside the frozen lists: {sorted(unlisted)}")

    output.mkdir(mode=0o750, parents=True)
    columns = [
        "donor_id",
        "well_id",
        "steatosis_observed",
        "steatosis_numeric",
        "steatosis_categorical",
        "steatosis_status",
        "fibrosis_observed",
        "fibrosis_any",
        "fibrosis_categorical",
        "fibrosis_status",
    ]
    for field in COVARIATE_FIELDS:
        columns.extend([f"{field}_observed", field])
    columns.append("is_adult")
    _write_tsv(output / "donor_endpoints.tsv", columns, endpoints)
    _write_tsv(
        output / "endpoint_masks.tsv",
        ["donor_id", "steatosis_observed", "fibrosis_observed", "is_adult"],
        endpoints,
    )
    write_json_exclusive(output / "donor_endpoints.json", endpoints)
    receipt = {
        "dataset_id": "gse296875",
        "role": spec["role"],
        "source_sheet": spec["source_sheet"],
        "source_supplement_sha256": spec["source_supplement_sha256"],
        "sheet_dump_sha256": sha256_file(cells),
        "join_lock_artifacts_sha256": sha256_file(join_lock / "ARTIFACTS.json"),
        "steatosis_endpoint": spec["steatosis"]["endpoint"],
        "steatosis_metric": spec["steatosis"]["metric"],
        "fibrosis_endpoint": spec["fibrosis"]["endpoint"],
        "fibrosis_metric": spec["fibrosis"]["metric"],
        "fibrosis_negative_source_values": list(spec["fibrosis"]["negative_source_values"]),
        "fibrosis_positive_source_values": list(spec["fibrosis"]["positive_source_values"]),
        "missingness": spec["inference"]["missingness"],
        "missing_encoding": "typed_null_with_observed_flag_never_zero_or_stage_zero",
        "forbidden_claims": list(spec["claims"]["forbidden"]),
        **summary,
    }
    write_json_exclusive(output / "endpoint_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse296875_phenotype_endpoint_lock",
            "dataset_id": "gse296875",
            "role": spec["role"],
            "unit_of_inference": "donor",
            "analyzed_donors": summary["analyzed_donors"],
            "steatosis_observed_donors": summary["steatosis_observed_donors"],
            "fibrosis_observed_donors": summary["fibrosis_observed_donors"],
            "missingness": spec["inference"]["missingness"],
            "ordinal_stage_inference": False,
            "recoding": "none",
            "champion_eligible": False,
            "external_or_sealed": False,
        },
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--join-lock", required=True, type=Path)
    parser.add_argument("--taskspec", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run(
        cells=arguments.cells,
        join_lock=arguments.join_lock,
        taskspec=arguments.taskspec,
        output=arguments.output,
    )
    print(json.dumps({"output": arguments.output.as_posix(), "receipt": receipt}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
