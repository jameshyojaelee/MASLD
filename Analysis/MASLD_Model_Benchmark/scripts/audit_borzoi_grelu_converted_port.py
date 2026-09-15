#!/usr/bin/env python3
"""Statically audit exact Borzoi gReLU converted-port PyTorch archives."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import pickletools
from typing import Any, Mapping
import zipfile


class ConvertedPortAuditError(RuntimeError):
    """Raised when converted-port sources differ from the frozen requirements."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def validate_file(path: Path, expected: Mapping[str, Any]) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ConvertedPortAuditError(f"source is not a regular file: {path.name}")
    observed = {
        "filename": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": digest(path),
    }
    if (
        observed["size_bytes"] != expected["size_bytes"]
        or observed["sha256"] != expected["sha256"]
    ):
        raise ConvertedPortAuditError(f"source bytes differ: {path.name}")
    return observed


def scan_pickle(
    payload: bytes,
    *,
    allowed_globals: set[str],
    forbidden_opcodes: set[str],
) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    globals_seen: set[str] = set()
    try:
        operations = pickletools.genops(payload)
        for operation, argument, _position in operations:
            counts[operation.name] += 1
            if operation.name in forbidden_opcodes:
                raise ConvertedPortAuditError(
                    f"forbidden pickle opcode: {operation.name}"
                )
            if operation.name == "GLOBAL":
                globals_seen.add(str(argument))
    except ValueError as error:
        raise ConvertedPortAuditError("pickle opcode stream is malformed") from error
    unexpected = globals_seen - allowed_globals
    if unexpected:
        raise ConvertedPortAuditError(
            f"unexpected pickle globals: {sorted(unexpected)}"
        )
    return {
        "pickle_opcode_counts": dict(sorted(counts.items())),
        "pickle_globals": sorted(globals_seen),
        "pickle_payload_unpickled": False,
    }


def audit_archive(
    path: Path,
    expected: Mapping[str, Any],
    archive_contract: Mapping[str, Any],
) -> dict[str, Any]:
    file_receipt = validate_file(path, expected)
    with path.open("rb") as handle:
        magic = handle.read(8)
    if not magic.startswith(b"PK\x03\x04") or magic.startswith(b"\x89HDF\r\n\x1a\n"):
        raise ConvertedPortAuditError("checkpoint archive magic differs")
    try:
        with zipfile.ZipFile(path, mode="r") as archive:
            infos = archive.infolist()
            if len(infos) > archive_contract["maximum_zip_entries_per_member"]:
                raise ConvertedPortAuditError("ZIP entry ceiling exceeded")
            names = [info.filename for info in infos]
            if len(names) != len(set(names)):
                raise ConvertedPortAuditError("duplicate ZIP member name")
            prefix = expected["archive_prefix"]
            for info in infos:
                item = PurePosixPath(info.filename)
                if (
                    item.is_absolute()
                    or ".." in item.parts
                    or not info.filename.startswith(f"{prefix}/")
                    or info.is_dir()
                    or info.flag_bits & 0x1
                    or info.compress_type != zipfile.ZIP_STORED
                ):
                    raise ConvertedPortAuditError(
                        f"unsafe or unexpected ZIP member: {info.filename}"
                    )
            pickle_names = [name for name in names if name == f"{prefix}/data.pkl"]
            storage_names = [
                name
                for name in names
                if name.startswith(f"{prefix}/data/")
                and name.removeprefix(f"{prefix}/data/").isdigit()
            ]
            if (
                len(infos) != archive_contract["zip_entries_per_member"]
                or len(pickle_names) != 1
                or len(storage_names)
                != archive_contract["tensor_storage_entries_per_member"]
                or sorted(int(name.rsplit("/", 1)[1]) for name in storage_names)
                != list(range(len(storage_names)))
                or f"{prefix}/byteorder" not in names
                or f"{prefix}/version" not in names
                or f"{prefix}/.data/serialization_id" not in names
            ):
                raise ConvertedPortAuditError("PyTorch ZIP member topology differs")
            data_pickle = archive.read(pickle_names[0])
            if (
                len(data_pickle) != archive_contract["data_pickle_size_bytes"]
                or sha256(data_pickle).hexdigest()
                != archive_contract["data_pickle_sha256"]
                or not data_pickle.startswith(b"\x80\x02")
            ):
                raise ConvertedPortAuditError("data.pkl identity or protocol differs")
            pickle_receipt = scan_pickle(
                data_pickle,
                allowed_globals=set(archive_contract["allowed_pickle_globals"]),
                forbidden_opcodes=set(archive_contract["forbidden_pickle_opcodes"]),
            )
            for key, value in archive_contract["pickle_opcode_counts"].items():
                if pickle_receipt["pickle_opcode_counts"].get(key) != value:
                    raise ConvertedPortAuditError(
                        f"pickle opcode count differs: {key}"
                    )
            if pickle_receipt["pickle_globals"] != sorted(
                archive_contract["allowed_pickle_globals"]
            ):
                raise ConvertedPortAuditError("pickle global roster differs")
            bad_member = archive.testzip()
            if bad_member is not None:
                raise ConvertedPortAuditError(
                    f"ZIP CRC failed for member: {bad_member}"
                )
            byteorder = archive.read(f"{prefix}/byteorder").decode("ascii")
            version = archive.read(f"{prefix}/version").decode("ascii").strip()
    except (OSError, zipfile.BadZipFile) as error:
        raise ConvertedPortAuditError("checkpoint is not a valid ZIP archive") from error
    file_receipt.update(
        {
            "replicate": expected["replicate"],
            "observed_file_magic_hex": magic[:4].hex(),
            "hdf5_signature_present": False,
            "archive_format": "pytorch_zip_serialization_with_pickle_metadata",
            "archive_prefix": expected["archive_prefix"],
            "zip_entry_count": len(infos),
            "tensor_storage_entry_count": len(storage_names),
            "data_pickle_size_bytes": len(data_pickle),
            "data_pickle_sha256": sha256(data_pickle).hexdigest(),
            "byteorder": byteorder,
            "serialization_version": version,
            "all_zip_member_crc_checks_passed": True,
            "archive_members_extracted": False,
            "tensor_storage_values_interpreted": False,
            **pickle_receipt,
        }
    )
    return file_receipt


def validate_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("schema_version") != (
        "masld-bench-borzoi-grelu-converted-port-static-audit-v2"
    ):
        raise ConvertedPortAuditError("contract schema differs")
    identity = contract["model_identity"]
    incident = contract["source_format_incident"]
    disposition = contract["admission_disposition"]
    safety = contract["safety_contract"]
    if (
        identity["model_id"] != "borzoi_grelu_converted_port"
        or identity["native_borzoi_substitute"] is not False
        or incident["observed_format"]
        != "PyTorch_ZIP_serialization_with_data.pkl_protocol_2"
        or incident["hdf5_signature_present"] is not False
        or incident["notebook_loader"] != "torch.load"
        or disposition["static_zip_and_pickle_opcode_scan_allowed"] is not True
        or disposition["unrestricted_torch_load_allowed"] is not False
        or disposition["torch_weights_only_load_allowed_by_this_contract"] is not False
        or disposition["model_forward_allowed"] is not False
        or disposition["global_frozen_census_mutation_allowed"] is not False
        or disposition["open_champion_eligible"] is not False
        or disposition["external_champion_claim_allowed"] is not False
        or disposition["universal_claim_allowed"] is not False
        or safety["pickle_payload_unpickled"] is not False
        or safety["tensor_storage_values_interpreted"] is not False
        or safety["archive_members_extracted"] is not False
        or safety["sealed_sources_accessed"] is not False
        or safety["biological_outcomes_accessed"] is not False
    ):
        raise ConvertedPortAuditError("fail-closed static-audit contract opened")
    members = contract["checkpoint_members"]
    if [member["replicate"] for member in members] != [0, 1, 2, 3]:
        raise ConvertedPortAuditError("replicate order differs")
    if [member["archive_prefix"] for member in members] != [
        "fold0",
        "fold1",
        "fold2",
        "fold3",
    ]:
        raise ConvertedPortAuditError("archive prefixes differ")
    if sum(member["size_bytes"] for member in members) != safety[
        "expected_total_download_bytes"
    ]:
        raise ConvertedPortAuditError("download byte total differs")


def audit(contract_path: Path, source_dir: Path) -> dict[str, Any]:
    if contract_path.is_symlink() or source_dir.is_symlink():
        raise ConvertedPortAuditError("symlinked contract or source directory prohibited")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    validate_contract(contract)
    small_sources = {
        name: validate_file(source_dir / name, expected)
        for name, expected in contract["small_sources"].items()
    }
    readme = (source_dir / "README.md").read_text(encoding="utf-8")
    notebook = (source_dir / "save_wandb_borzoi_ckpt_human.ipynb").read_text(
        encoding="utf-8"
    )
    if (
        "license: mit" not in readme
        or "HDF5" not in readme
        or "human_state_dict_rep0.h5" not in readme
        or 'torch.load(f\\"/data/borzoi/torch_weights/fold{rep}.h5\\")'
        not in notebook
        or "lm.model.load_state_dict(state_dict)" not in notebook
    ):
        raise ConvertedPortAuditError("model-card or conversion-notebook semantics differ")
    members = [
        audit_archive(
            source_dir / expected["filename"],
            expected,
            {
                **contract["static_archive_expectations"],
                **contract["safety_contract"],
            },
        )
        for expected in contract["checkpoint_members"]
    ]
    if len({member["data_pickle_sha256"] for member in members}) != 1:
        raise ConvertedPortAuditError("replicate state-dictionary schemas differ")
    return {
        "schema_version": "masld-bench-borzoi-grelu-converted-port-audit-v2",
        "model_id": contract["model_identity"]["model_id"],
        "source_revision": contract["model_identity"]["source_revision"],
        "small_sources": small_sources,
        "source_format_incident": contract["source_format_incident"],
        "checkpoint_members": members,
        "checkpoint_bytes_present": True,
        "pickle_payload_unpickled": False,
        "unrestricted_torch_load_executed": False,
        "torch_weights_only_load_executed": False,
        "model_forward_executed": False,
        "native_borzoi_substitute": False,
        "open_champion_eligible": False,
        "external_champion_claim_allowed": False,
        "universal_claim_allowed": False,
        "terminal_disposition": (
            "static_archive_scan_passed_if_frozen_but_separate_weights_only_load_"
            "runtime_native_parity_and_weight_authority_gates_remain_closed"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ConvertedPortAuditError("refusing to overwrite audit output")
    receipt = audit(args.contract, args.source_dir)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "converted_port_static_audit.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
