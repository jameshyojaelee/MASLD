#!/usr/bin/env python3
"""Outcome-blind exact verification of the Yakubovsky common gene axis.

This command opens only the top-level donor references, patient strings, and
``gene_name`` character objects.  It never reads ``mat`` or any lipid field.
The static certificate is safe to reuse only because it is bound to the exact
SHA256 of ``v.mat``.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from yakubovsky_common import (
    ContractError,
    RELEASE_ID,
    atomic_write_frame,
    bulk_decode_matlab_cellstr,
    decode_matlab_char,
    default_paths,
    gene_axis_sha256,
    sha256_file,
    validate_gene_axis_manifest,
)


def verify(paths: dict[str, Path], manifest_path: Path) -> dict[str, object]:
    v_mat_sha256 = sha256_file(paths["v_mat"])
    contract = validate_gene_axis_manifest(manifest_path, v_mat_sha256)

    with h5py.File(paths["v_mat"], "r") as handle:
        references = np.asarray(handle["v"][...]).reshape(-1, order="F")
        if len(references) != len(contract.rows):
            raise ContractError(
                f"Gene-axis donor count drift: {len(references)} != {len(contract.rows)}"
            )
        axes: dict[str, h5py.Dataset] = {}
        for reference, expected in zip(
            references, contract.rows.itertuples(index=False), strict=True
        ):
            group = handle[reference]
            donor = decode_matlab_char(group["patient"])
            if donor != expected.donor or group.name != expected.group_path:
                raise ContractError(
                    "Gene-axis donor/group order drift: "
                    f"{donor}@{group.name} != {expected.donor}@{expected.group_path}"
                )
            axes[donor] = group["gene_name"]
        resolved = bulk_decode_matlab_cellstr(axes, handle, allow_null=False)

    observed_hashes: dict[str, str] = {}
    for row in contract.rows.itertuples(index=False):
        genes = resolved.values[row.donor]
        if len(genes) != int(row.n_genes):
            raise ContractError(
                f"Gene-axis length drift for {row.donor}: {len(genes)} != {row.n_genes}"
            )
        observed = gene_axis_sha256(genes)
        if observed != row.gene_axis_sha256:
            raise ContractError(
                f"Gene-axis content drift for {row.donor}: "
                f"{observed} != {row.gene_axis_sha256}"
            )
        observed_hashes[row.donor] = observed
    if set(observed_hashes.values()) != {contract.common_sha256}:
        raise ContractError("Exact all-donor gene-axis equality did not rederive")

    producer = Path(__file__).resolve()
    return {
        "release_id": RELEASE_ID,
        "status": "exact_all_donor_gene_axes_verified",
        "n_donors": len(observed_hashes),
        "n_genes": contract.n_genes,
        "common_gene_axis_sha256": contract.common_sha256,
        "v_mat_sha256": contract.v_mat_sha256,
        "manifest_sha256": contract.manifest_sha256,
        "bulk_reference_count": resolved.n_references,
        "bulk_unique_character_objects": resolved.n_unique_objects,
        "bulk_objects_visited": resolved.n_objects_visited,
        "producer": str(producer),
        "producer_sha256": sha256_file(producer),
        "python": sys.version.replace("\n", " "),
        "verified_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "matrix_values_read": False,
        "lipid_fields_read": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    paths = default_paths(args.base)
    manifest_path = (args.manifest or paths["gene_axis_manifest"]).resolve()
    output = (args.output or paths["gene_axis_verification"]).resolve()
    result = verify(paths, manifest_path)
    stable_fields = [
        "release_id",
        "status",
        "n_donors",
        "n_genes",
        "common_gene_axis_sha256",
        "v_mat_sha256",
        "manifest_sha256",
        "bulk_reference_count",
        "bulk_unique_character_objects",
        "producer",
        "producer_sha256",
        "matrix_values_read",
        "lipid_fields_read",
    ]
    if args.write:
        if output.exists():
            raise ContractError(f"Refusing to overwrite gene-axis verification: {output}")
        atomic_write_frame(output, pd.DataFrame([result]))
        result["verification_artifact"] = str(output)
        result["verification_artifact_sha256"] = sha256_file(output)
    else:
        if not output.is_file():
            raise ContractError(f"Missing gene-axis verification artifact: {output}")
        stored = pd.read_csv(output, sep="\t", dtype=str, keep_default_na=False)
        if len(stored) != 1:
            raise ContractError("Gene-axis verification artifact must have exactly one row")
        for field in stable_fields:
            expected = str(result[field])
            observed = stored.iloc[0].get(field, "")
            if observed.lower() in {"true", "false"}:
                observed = observed.title()
            if observed != expected:
                raise ContractError(
                    f"Gene-axis verification artifact drift for {field}: "
                    f"{observed!r} != {expected!r}"
                )
        result["verification_artifact"] = str(output)
        result["verification_artifact_sha256"] = sha256_file(output)
        result["stored_artifact_rederived"] = True
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
