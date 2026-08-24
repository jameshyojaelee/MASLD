#!/usr/bin/env python3
"""Freeze the exposed, replicate-safe GSE281364 MPRA count substrate."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path


class OutcomeAuditError(ValueError):
    """Raised when the exposed MPRA replicate topology differs."""


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _expected_alleles(path: Path) -> dict[str, set[str]]:
    output: dict[str, set[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            construct = line.split("\t", 1)[0]
            if not construct:
                raise OutcomeAuditError("empty construct in barcode mapping")
            alternative = construct.endswith("_Mut")
            element_id = construct[:-4] if alternative else construct
            output.setdefault(element_id, set()).add("alt" if alternative else "ref")
    if (
        len(output) != 5_442
        or sum(alleles == {"ref", "alt"} for alleles in output.values()) != 5_355
        or sum(alleles == {"ref"} for alleles in output.values()) != 81
        or sum(alleles == {"alt"} for alleles in output.values()) != 6
    ):
        raise OutcomeAuditError("barcode-mapping allele topology differs")
    return output


def audit(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists() or any(
        path.is_symlink()
        for path in (
            arguments.hepg2_counts,
            arguments.lx2_counts,
            arguments.mapping_counts,
            arguments.gate,
        )
    ):
        raise OutcomeAuditError("MPRA outcome admission request differs")
    arguments.output.mkdir(mode=0o750)
    gate = json.loads(arguments.gate.read_text())
    if (
        not gate.get("source_reproduction_complete")
        or not gate.get("source_reproduction_pass")
        or gate.get("official_source_calls_authoritative") is not True
    ):
        raise OutcomeAuditError("MPRA source-reproduction gate differs")
    conditions = {
        "HepG2": {"control", "PAOA"},
        "LX2": {"control", "TGFb"},
    }
    expected_alleles = _expected_alleles(arguments.mapping_counts)
    observed_rows: dict[tuple[str, str, str], dict[str, object]] = {}
    keys: set[tuple[str, str, str]] = set()
    sample_context: dict[str, str] = {}
    element_alleles: dict[str, set[str]] = {}
    for cell_line, path in (
        ("HepG2", arguments.hepg2_counts),
        ("LX2", arguments.lx2_counts),
    ):
        rows = _read_rows(path)
        for row in rows:
            condition = row["condition"]
            replicate = int(row["replicate"])
            allele = row["allele"]
            sample_id = row["sample_id"]
            element_id = row["element_id"]
            if (
                condition not in conditions[cell_line]
                or replicate not in {1, 2, 3, 4}
                or allele not in {"ref", "alt"}
                or element_id not in expected_alleles
                or allele not in expected_alleles[element_id]
            ):
                raise OutcomeAuditError("MPRA replicate row differs")
            context_id = f"{cell_line}_{condition}"
            previous_context = sample_context.setdefault(sample_id, context_id)
            if previous_context != context_id:
                raise OutcomeAuditError("sample belongs to multiple MPRA contexts")
            key = (element_id, allele, sample_id)
            if key in keys:
                raise OutcomeAuditError("duplicate MPRA replicate row")
            keys.add(key)
            element_alleles.setdefault(element_id, set()).add(allele)
            dna, rna, barcodes = (int(row[field]) for field in ("DNA", "RNA", "n_barcodes"))
            if min(dna, rna, barcodes) < 0:
                raise OutcomeAuditError("negative MPRA count")
            observed_rows[key] = {
                "element_id": element_id,
                "allele": allele,
                "context_id": context_id,
                "cell_line": cell_line,
                "condition": condition,
                "experimental_replicate": replicate,
                "sample_id": sample_id,
                "DNA": dna,
                "RNA": rna,
                "n_barcodes": barcodes,
                "assay_state": "observed",
                "missing_reason": "not_applicable",
                "pairing": "same_sample_different_aliquot",
                "biological_unit": "experimental_replicate",
                "donor_id": "not_applicable",
            }
    context_samples: dict[str, set[str]] = {}
    for sample, context in sample_context.items():
        context_samples.setdefault(context, set()).add(sample)
    if set(context_samples) != {
        "HepG2_control",
        "HepG2_PAOA",
        "LX2_control",
        "LX2_TGFb",
    } or any(len(samples) != 4 for samples in context_samples.values()):
        raise OutcomeAuditError("MPRA context/sample topology differs")
    output_rows: list[dict[str, object]] = []
    missing_rows = 0
    for sample_id, context_id in sorted(sample_context.items()):
        cell_line, condition = context_id.split("_", 1)
        replicate = int(sample_id.rsplit("_r", 1)[1])
        for element_id, alleles in sorted(expected_alleles.items()):
            for allele in sorted(alleles):
                key = (element_id, allele, sample_id)
                observed = observed_rows.get(key)
                if observed is not None:
                    output_rows.append(observed)
                    continue
                missing_rows += 1
                output_rows.append(
                    {
                        "element_id": element_id,
                        "allele": allele,
                        "context_id": context_id,
                        "cell_line": cell_line,
                        "condition": condition,
                        "experimental_replicate": replicate,
                        "sample_id": sample_id,
                        "DNA": "",
                        "RNA": "",
                        "n_barcodes": "",
                        "assay_state": "below_qc",
                        "missing_reason": "absent_from_source_reproduction_aggregate",
                        "pairing": "same_sample_different_aliquot",
                        "biological_unit": "experimental_replicate",
                        "donor_id": "not_applicable",
                    }
                )
    if len(observed_rows) != 172_528 or len(output_rows) != 172_752 or missing_rows != 224:
        raise OutcomeAuditError(
            "combined observed/explicit-missing MPRA census differs"
        )
    output_rows.sort(
        key=lambda row: (
            str(row["element_id"]),
            str(row["allele"]),
            str(row["context_id"]),
            int(row["experimental_replicate"]),
        )
    )
    table_path = arguments.output / "replicate_outcomes.tsv.gz"
    with table_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(
                    text,
                    fieldnames=list(output_rows[0]),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(output_rows)
    receipt = {
        "schema_version": "masld-bench-gse281364-replicate-outcomes-v1",
        "status": "pass",
        "rows": len(output_rows),
        "elements": len(expected_alleles),
        "paired_construct_elements": sum(
            alleles == {"ref", "alt"} for alleles in expected_alleles.values()
        ),
        "reference_only_construct_elements": sum(
            alleles == {"ref"} for alleles in expected_alleles.values()
        ),
        "alternative_only_construct_elements": sum(
            alleles == {"alt"} for alleles in expected_alleles.values()
        ),
        "observed_rows": len(observed_rows),
        "below_qc_rows": missing_rows,
        "missing_evidence_encoded_as_zero": False,
        "contexts": {key: sorted(value) for key, value in context_samples.items()},
        "experimental_replicates": len(sample_context),
        "donor_count": 0,
        "replicates_are_independent_donors": False,
        "pairing": "same_sample_different_aliquot_DNA_RNA",
        "activity_transform_fit_scope": "not_fit_source_admission_only",
        "outcome_role": "exposed_development_MPRA_only",
        "eQTL_or_ieQTL_supervision": False,
        "sealed_outcomes_loaded": False,
        "champion_eligible": False,
        "replicate_outcomes_sha256": _sha256_file(table_path),
    }
    (arguments.output / "outcome_admission_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hepg2-counts", type=Path, required=True)
    parser.add_argument("--lx2-counts", type=Path, required=True)
    parser.add_argument("--mapping-counts", type=Path, required=True)
    parser.add_argument("--gate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    audit(parser.parse_args())


if __name__ == "__main__":
    main()
