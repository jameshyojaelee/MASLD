#!/usr/bin/env python3
"""Freeze and hash the well-namespaced GSE296875 barcode-to-donor join.

The authoritative join lives outside the benchmark tree and is read-only there.
This campaign copies it in under a checksum, proves the properties the dataset
requirements depend on, and records the ones the requirements never declared: whether
a bare 10x barcode recurs across wells, and whether a donor spans wells.  A
nucleus that cannot be resolved to a donor is absent from the roster, never
assigned a default donor.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import shutil
from typing import Any

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


REQUIRED_COLUMNS = ("cell_id", "donor_id", "well_id", "raw_barcode")
WELL_PREFIX = re.compile(r"^(well[1-9][0-9]*)_(.+)$")


class BarcodeJoinError(ValueError):
    """Raised when the barcode-to-donor join contradicts the frozen requirements."""


def read_join(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or tuple(reader.fieldnames) != REQUIRED_COLUMNS:
            raise BarcodeJoinError(f"join columns differ: {reader.fieldnames}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise BarcodeJoinError("join is empty")
    return rows


def audit_join(rows: list[dict[str, str]]) -> dict[str, Any]:
    """Derive every join property without recoding or filling a missing donor."""

    cell_ids: set[str] = set()
    namespaced: set[str] = set()
    donor_wells: dict[str, set[str]] = {}
    well_donors: dict[str, set[str]] = {}
    well_nuclei: dict[str, int] = {}
    bare_wells: dict[str, set[str]] = {}
    bare_counts: dict[str, int] = {}
    blank_donor_rows = 0
    for row in rows:
        cell_id = row["cell_id"].strip()
        donor_id = row["donor_id"].strip()
        well_id = row["well_id"].strip()
        raw_barcode = row["raw_barcode"].strip()
        if not cell_id or not well_id or not raw_barcode:
            raise BarcodeJoinError(f"join row has a blank key field: {row}")
        if not donor_id:
            blank_donor_rows += 1
            continue
        if cell_id in cell_ids:
            raise BarcodeJoinError(f"cell_id is not unique: {cell_id}")
        cell_ids.add(cell_id)
        match = WELL_PREFIX.fullmatch(raw_barcode)
        if match is None:
            raise BarcodeJoinError(f"raw_barcode is not well-namespaced: {raw_barcode}")
        prefix, bare = match.group(1), match.group(2)
        if prefix != well_id:
            raise BarcodeJoinError(
                f"raw_barcode well prefix disagrees with well_id: {raw_barcode} vs {well_id}"
            )
        if cell_id != f"{well_id}_{bare}":
            raise BarcodeJoinError(
                f"cell_id is not the well-namespaced barcode: {cell_id} vs {well_id}_{bare}"
            )
        namespaced.add(raw_barcode)
        donor_wells.setdefault(donor_id, set()).add(well_id)
        well_donors.setdefault(well_id, set()).add(donor_id)
        well_nuclei[well_id] = well_nuclei.get(well_id, 0) + 1
        bare_wells.setdefault(bare, set()).add(well_id)
        bare_counts[bare] = bare_counts.get(bare, 0) + 1
    if blank_donor_rows:
        raise BarcodeJoinError(
            f"{blank_donor_rows} nuclei carry no donor; a join row may not default a donor"
        )
    if len(namespaced) != len(cell_ids):
        raise BarcodeJoinError("well-namespaced barcodes are not unique")
    recurring = {bare: wells for bare, wells in bare_wells.items() if len(wells) > 1}
    multi_well_donors = {
        donor: sorted(wells) for donor, wells in donor_wells.items() if len(wells) > 1
    }
    return {
        "nuclei": len(cell_ids),
        "donors": len(donor_wells),
        "wells": len(well_nuclei),
        "unique_cell_id": len(cell_ids),
        "unique_well_namespaced_barcode": len(namespaced),
        "unique_bare_10x_barcode": len(bare_counts),
        "bare_10x_barcodes_shared_across_wells": len(recurring),
        "maximum_wells_sharing_one_bare_10x_barcode": (
            max((len(wells) for wells in bare_wells.values()), default=0)
        ),
        "nuclei_whose_bare_barcode_is_not_globally_unique": sum(
            bare_counts[bare] for bare in recurring
        ),
        "bare_barcode_is_a_unique_key": not recurring,
        "well_namespacing_is_load_bearing": bool(recurring),
        "donors_spanning_more_than_one_well": len(multi_well_donors),
        "multi_well_donors": multi_well_donors,
        "donor_is_nested_within_well": not multi_well_donors,
        "well_nuclei": dict(sorted(well_nuclei.items())),
        "well_donors": {well: sorted(donors) for well, donors in sorted(well_donors.items())},
        "donor_well": {donor: sorted(wells)[0] for donor, wells in sorted(donor_wells.items())},
    }


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    payload = "\t".join(header) + "\n"
    payload += "".join("\t".join(row) + "\n" for row in rows)
    path.write_text(payload, encoding="utf-8")
    path.chmod(0o440)


def run(
    *,
    join: Path,
    source_lock: Path,
    raw_rna_sources: Path,
    membership: Path,
    expected_donors: int,
    expected_wells: int,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise BarcodeJoinError("refusing to overwrite barcode-join artifact")
    source_digest = sha256_file(join)
    rows = read_join(join)
    audit = audit_join(rows)
    if audit["donors"] != expected_donors:
        raise BarcodeJoinError(
            f"donor count differs: observed {audit['donors']}, frozen {expected_donors}"
        )
    if audit["wells"] != expected_wells:
        raise BarcodeJoinError(
            f"well count differs: observed {audit['wells']}, frozen {expected_wells}"
        )

    lock = json.loads(source_lock.read_text(encoding="utf-8"))
    declared_wells = sorted(lock["raw_rna_sources"])
    if declared_wells != sorted(audit["well_nuclei"]):
        raise BarcodeJoinError(
            f"source-lock wells differ from the join: {declared_wells} vs "
            f"{sorted(audit['well_nuclei'])}"
        )
    with raw_rna_sources.open(encoding="utf-8", newline="") as handle:
        rna_rows = list(csv.DictReader(handle, delimiter="\t"))
    rna_gsm = {row["well"]: row["gsm"] for row in rna_rows}
    if rna_gsm != {well: lock["raw_rna_sources"][well]["gsm"] for well in declared_wells}:
        raise BarcodeJoinError("raw RNA GSM manifest disagrees with the source lock")

    membership_contract = json.loads((membership / "contract.json").read_text(encoding="utf-8"))
    if membership_contract["mapping_sha256"] != source_digest:
        raise BarcodeJoinError(
            "the frozen fragment-membership campaign used a different join checksum"
        )
    with (membership / "well_census.tsv").open(encoding="utf-8", newline="") as handle:
        membership_wells = list(csv.DictReader(handle, delimiter="\t"))
    for row in membership_wells:
        well = row["well_id"]
        if int(row["nuclei"]) != audit["well_nuclei"][well]:
            raise BarcodeJoinError(f"nuclei per well disagree with the frozen census: {well}")
        if int(row["donors"]) != len(audit["well_donors"][well]):
            raise BarcodeJoinError(f"donors per well disagree with the frozen census: {well}")
    atac_gsm = {row["well_id"]: row["gsm"] for row in membership_wells}

    output.mkdir(mode=0o750, parents=True)
    copied = output / join.name
    shutil.copyfile(join, copied)
    copied.chmod(0o440)
    if sha256_file(copied) != source_digest:
        raise BarcodeJoinError("the in-tree join copy is not byte-identical to its source")
    _write_tsv(
        output / "well_census.tsv",
        ["well_id", "rna_gsm", "atac_gsm", "nuclei", "donors"],
        [
            [
                well,
                rna_gsm[well],
                atac_gsm.get(well, ""),
                str(audit["well_nuclei"][well]),
                str(len(audit["well_donors"][well])),
            ]
            for well in sorted(audit["well_nuclei"])
        ],
    )
    _write_tsv(
        output / "donor_census.tsv",
        ["donor_id", "well_id", "nuclei"],
        [
            [donor, well, str(sum(1 for row in rows if row["donor_id"].strip() == donor))]
            for donor, well in audit["donor_well"].items()
        ],
    )
    receipt = {
        "dataset_id": "gse296875",
        "biological_unit": "donor",
        "technical_batch": "well",
        "join_source_path": join.as_posix(),
        "join_sha256": source_digest,
        "join_size_bytes": join.stat().st_size,
        "source_lock_sha256": sha256_file(source_lock),
        "raw_rna_sources_sha256": sha256_file(raw_rna_sources),
        "fragment_membership_artifacts_sha256": sha256_file(membership / "ARTIFACTS.json"),
        "expected_donors": expected_donors,
        "expected_wells": expected_wells,
        "missing_barcode_encoding": "absent_from_roster_never_a_default_donor",
        **{key: value for key, value in audit.items() if key not in {"well_donors"}},
    }
    write_json_exclusive(output / "join_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse296875_barcode_donor_join_lock",
            "dataset_id": "gse296875",
            "biological_unit": "donor",
            "technical_batch": "well",
            "join_sha256": source_digest,
            "donors": audit["donors"],
            "wells": audit["wells"],
            "nuclei": audit["nuclei"],
            "donor_is_nested_within_well": audit["donor_is_nested_within_well"],
            "well_namespacing_is_load_bearing": audit["well_namespacing_is_load_bearing"],
            "champion_eligible": False,
            "external_or_sealed": False,
        },
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--join", required=True, type=Path)
    parser.add_argument("--source-lock", required=True, type=Path)
    parser.add_argument("--raw-rna-sources", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--expected-donors", required=True, type=int)
    parser.add_argument("--expected-wells", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run(
        join=arguments.join,
        source_lock=arguments.source_lock,
        raw_rna_sources=arguments.raw_rna_sources,
        membership=arguments.membership,
        expected_donors=arguments.expected_donors,
        expected_wells=arguments.expected_wells,
        output=arguments.output,
    )
    print(json.dumps({"output": arguments.output.as_posix(), "receipt": receipt}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
