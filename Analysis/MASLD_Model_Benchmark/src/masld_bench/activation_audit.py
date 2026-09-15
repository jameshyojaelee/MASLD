"""Read-only audits for admitting existing scientific data authorities."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping


class ActivationAuditError(ValueError):
    """Raised when an upstream dataset authority cannot be verified."""


_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_ENSEMBL_GENE = re.compile(r"^ENSG[0-9]{11}(?:\.[0-9]+)?$")


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _write_json_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    payload = json.dumps(
        value,
        sort_keys=True,
        indent=2,
        ensure_ascii=False,
        allow_nan=False,
    ) + "\n"
    with path.open("x", encoding="utf-8") as handle:
        handle.write(payload)
    path.chmod(0o440)


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ActivationAuditError(f"{label} must be a JSON object")
    return value


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ActivationAuditError(f"cannot read contract lock: {error}") from error
    contract = dict(_mapping(raw, "contract lock"))
    required = {
        "schema_version",
        "atlas_realpath",
        "atlas_size_bytes",
        "cell_order_sha256",
        "config_sha256",
        "gene_order_sha256",
        "lock_sha256",
        "manifest_sha256",
        "observed_contract",
        "production_ready",
        "raw_counts",
        "raw_counts_sha256",
    }
    if set(contract) != required:
        raise ActivationAuditError("contract lock has an unexpected schema")
    if contract["schema_version"] != "masld-cl-contract-v1":
        raise ActivationAuditError("contract lock schema_version is not supported")
    if contract["production_ready"] is not True:
        raise ActivationAuditError("upstream atlas contract is not production-ready")
    for field in (
        "cell_order_sha256",
        "config_sha256",
        "gene_order_sha256",
        "lock_sha256",
        "raw_counts_sha256",
    ):
        if not isinstance(contract[field], str) or not _SHA256.fullmatch(contract[field]):
            raise ActivationAuditError(f"contract {field} is not a SHA-256 digest")
    identity = dict(contract)
    claimed = identity.pop("lock_sha256")
    if _canonical_sha256(identity) != claimed:
        raise ActivationAuditError("contract lock identity hash does not rederive")
    return contract


def _decode_attribute(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, list):
        return [_decode_attribute(item) for item in value]
    return value


def _dataset_length(group: Any, key: str) -> int | None:
    if key not in group:
        return None
    value = group[key]
    if not hasattr(value, "shape") or len(value.shape) < 1:
        return None
    return int(value.shape[0])


def _string_values(dataset: Any) -> list[str]:
    try:
        values = dataset.asstr()[:]
    except (AttributeError, TypeError, ValueError):
        values = dataset[:]
    result: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            result.append(value.decode("utf-8"))
        else:
            result.append(str(value))
    return result


def _categorical_levels(group: Any, key: str) -> list[str] | None:
    if key not in group:
        return None
    value = group[key]
    if hasattr(value, "keys") and "categories" in value:
        return _string_values(value["categories"])
    return None


def audit_resource_atlas(
    *,
    atlas_path: str | Path,
    contract_lock_path: str | Path,
    output_dir: str | Path,
) -> Path:
    """Verify the upstream atlas requirements and inventory Geneformer inputs."""

    try:
        import h5py
    except ImportError as error:
        raise ActivationAuditError("h5py is required for the atlas audit") from error

    atlas = Path(atlas_path).resolve(strict=True)
    contract_path = Path(contract_lock_path).resolve(strict=True)
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    contract = _load_contract(contract_path)
    if Path(str(contract["atlas_realpath"])).resolve(strict=True) != atlas:
        raise ActivationAuditError("contract atlas path differs from the requested atlas")
    if atlas.stat().st_size != int(contract["atlas_size_bytes"]):
        raise ActivationAuditError("atlas size differs from the upstream contract")

    manifest_hashes = dict(_mapping(contract["manifest_sha256"], "manifest_sha256"))
    manifest_observations: dict[str, dict[str, Any]] = {}
    for filename, expected in sorted(manifest_hashes.items()):
        if not isinstance(expected, str) or not _SHA256.fullmatch(expected):
            raise ActivationAuditError(f"invalid manifest hash for {filename}")
        source = contract_path.parent / filename
        observed = _sha256_file(source)
        if observed != expected:
            raise ActivationAuditError(f"upstream manifest changed: {filename}")
        manifest_observations[filename] = {
            "path": source.as_posix(),
            "sha256": observed,
            "size_bytes": source.stat().st_size,
        }

    with h5py.File(atlas, "r") as handle:
        for key in ("obs", "raw", "raw/X", "raw/var"):
            if key not in handle:
                raise ActivationAuditError(f"atlas lacks required H5AD object {key!r}")
        obs = handle["obs"]
        raw = handle["raw"]
        raw_x = handle["raw/X"]
        raw_var = handle["raw/var"]
        obs_keys = sorted(map(str, obs.keys()))
        var_keys = sorted(map(str, raw_var.keys()))
        n_cells = _dataset_length(obs, "_index")
        n_genes = _dataset_length(raw_var, "_index")
        if n_cells is None or n_genes is None:
            raise ActivationAuditError("atlas obs/raw-var indices are unavailable")
        expected_shape = [n_cells, n_genes]
        contract_counts = dict(_mapping(contract["raw_counts"], "raw_counts"))
        if list(contract_counts.get("shape", ())) != expected_shape:
            raise ActivationAuditError("H5AD axes differ from the raw-count contract")
        matrix_shape = _decode_attribute(raw_x.attrs.get("shape"))
        if matrix_shape is not None and list(matrix_shape) != expected_shape:
            raise ActivationAuditError("raw/X shape attribute differs from its axes")

        gene_ids = _string_values(raw_var["_index"])
        ensembl_index_fraction = sum(
            _ENSEMBL_GENE.fullmatch(value) is not None for value in gene_ids
        ) / len(gene_ids)
        source_gene_ids = (
            _string_values(raw_var["gene_ids"])
            if "gene_ids" in raw_var
            else []
        )
        versionless_source_gene_ids = [
            value.split(".", 1)[0] for value in source_gene_ids
        ]
        source_gene_id_contract = {
            "field": "gene_ids" if source_gene_ids else None,
            "n_values": len(source_gene_ids),
            "ensembl_fraction": (
                sum(_ENSEMBL_GENE.fullmatch(value) is not None for value in source_gene_ids)
                / len(source_gene_ids)
                if source_gene_ids
                else 0.0
            ),
            "empty_values": sum(not value for value in source_gene_ids),
            "versionless_duplicate_count": (
                len(versionless_source_gene_ids)
                - len(set(versionless_source_gene_ids))
            ),
        }
        ensembl_id_derivable = (
            source_gene_id_contract["n_values"] == n_genes
            and source_gene_id_contract["ensembl_fraction"] == 1.0
            and source_gene_id_contract["empty_values"] == 0
            and source_gene_id_contract["versionless_duplicate_count"] == 0
        )
        geneformer_fields = {
            "raw_unselected_counts": True,
            "obs_n_counts": "n_counts" in obs,
            "var_ensembl_id": "ensembl_id" in raw_var,
            "var_ensembl_id_derivable_from_gene_ids": ensembl_id_derivable,
        }
        likely_gene_identity_fields = [
            key
            for key in var_keys
            if any(token in key.lower() for token in ("ensembl", "gene", "symbol", "feature"))
        ]
        likely_label_fields = [
            key
            for key in obs_keys
            if any(token in key.lower() for token in ("cell_type", "celltype", "lineage", "label"))
        ]
        likely_group_fields = [
            key
            for key in obs_keys
            if any(token in key.lower() for token in ("donor", "sample", "library", "dataset"))
        ]
        categorical_levels = {
            key: levels
            for key in (
                "cell_type",
                "cell_type_raw",
                "dataset",
                "sample",
                "preparation_method",
            )
            if (levels := _categorical_levels(obs, key)) is not None
        }

    observed_contract = dict(_mapping(contract["observed_contract"], "observed_contract"))
    descriptive = dict(_mapping(observed_contract["descriptive"], "descriptive contract"))
    if n_cells != int(descriptive["cells"]):
        raise ActivationAuditError("atlas cell axis differs from the descriptive contract")
    tokenizer_ready = (
        geneformer_fields["raw_unselected_counts"]
        and geneformer_fields["obs_n_counts"]
        and (
            geneformer_fields["var_ensembl_id"]
            or geneformer_fields["var_ensembl_id_derivable_from_gene_ids"]
        )
    )
    blockers: list[str] = []
    if not geneformer_fields["obs_n_counts"]:
        blockers.append("derive_and_freeze_obs_n_counts_from_raw_X")
    if (
        not geneformer_fields["var_ensembl_id"]
        and geneformer_fields["var_ensembl_id_derivable_from_gene_ids"]
    ):
        blockers.append("freeze_versionless_gene_ids_to_ensembl_id_derivation")
    elif not geneformer_fields["var_ensembl_id"]:
        blockers.append("freeze_one_to_one_gencode_v49_ensembl_gene_crosswalk")
    if not likely_label_fields:
        blockers.append("freeze_cell_label_source_and_broad_liver_five_mapping")
    if not likely_group_fields:
        blockers.append("freeze_donor_library_and_cohort_group_fields")

    payload = {
        "schema_version": "masld-bench-resource-atlas-activation-audit-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_id": "resource_atlas_current",
        "source": {
            "atlas_path": atlas.as_posix(),
            "atlas_size_bytes": atlas.stat().st_size,
            "contract_lock_path": contract_path.as_posix(),
            "contract_lock_file_sha256": _sha256_file(contract_path),
            "contract_lock_id": contract["lock_sha256"],
            "raw_counts_sha256": contract["raw_counts_sha256"],
            "cell_order_sha256": contract["cell_order_sha256"],
            "gene_order_sha256": contract["gene_order_sha256"],
            "manifest_artifacts": manifest_observations,
        },
        "census": {
            "cells": n_cells,
            "genes": n_genes,
            "raw_counts_shape": expected_shape,
            "raw_counts_nnz": int(contract_counts["nnz"]),
            "descriptive_donors": int(descriptive["donors"]),
            "analyzed_donors": int(observed_contract["analyzed"]["donors"]),
        },
        "h5ad_inventory": {
            "obs_fields": obs_keys,
            "raw_var_fields": var_keys,
            "raw_var_index_ensembl_fraction": ensembl_index_fraction,
            "raw_var_gene_ids_contract": source_gene_id_contract,
            "likely_gene_identity_fields": likely_gene_identity_fields,
            "likely_label_fields": likely_label_fields,
            "likely_group_fields": likely_group_fields,
            "categorical_levels": categorical_levels,
        },
        "geneformer_input_contract": {
            **geneformer_fields,
            "ready_without_derivation": tokenizer_ready,
            "blockers": blockers,
        },
        "dataset_activation_ready": False,
        "dataset_activation_blockers": [
            *blockers,
            "freeze_rights_topology_join_labels_qc_reference_and_exposure_authorities",
            "freeze_donor_grouped_1000_cell_smoke_subset_without_sealed_outcomes",
        ],
    }
    payload["audit_sha256"] = _canonical_sha256(payload)
    destination = output / "resource_atlas_activation_audit.json"
    _write_json_exclusive(destination, payload)
    return destination


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    destination = audit_resource_atlas(
        atlas_path=args.atlas,
        contract_lock_path=args.contract_lock,
        output_dir=args.output,
    )
    print(destination.as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
