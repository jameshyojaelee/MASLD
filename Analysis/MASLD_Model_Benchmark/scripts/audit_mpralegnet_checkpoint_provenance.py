#!/usr/bin/env python3
"""Audit MPRALegNet author-checkpoint equivalence without reading outcomes."""

from __future__ import annotations

import argparse
import builtins
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
from hashlib import md5, sha256
import io
import json
import pickle
import pickletools
from pathlib import Path
import socket
import tomllib
from typing import Any, Mapping, Sequence
import zipfile

from safetensors.torch import load_file, save
import torch

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-mpralegnet-provenance-audit-v1"
STATUS = "prespecified_outcome_blind_provenance_equivalence_audit"
MODEL_ID = "mpralegnet_hepg2_test1_val2"
AUTHOR_MEMBER = (
    "final_dump/models/HepG2/md_shift_reverse_noavg_noch/"
    "best_model_test1_val2.ckpt"
)
AUTHOR_CONFIG_MEMBER = (
    "final_dump/models/HepG2/md_shift_reverse_noavg_noch/config.json"
)
HF_CHECKPOINT_SHA256 = (
    "470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39"
)


class ProvenanceAuditError(RuntimeError):
    """Raised when an authority or outcome-access separation differs."""


def file_digest(path: Path, algorithm: str = "sha256") -> str:
    value = sha256() if algorithm == "sha256" else md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def bytes_sha256(value: bytes) -> str:
    return sha256(value).hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ProvenanceAuditError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise ProvenanceAuditError(f"{label} must be an object")
    return value


def load_config(path: Path) -> dict[str, Any]:
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise ProvenanceAuditError(f"invalid audit config: {error}") from error
    if not isinstance(value, dict):
        raise ProvenanceAuditError("audit config must be a table")
    return value


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def safe_project_path(root: Path, relative: object, *, label: str) -> Path:
    if not isinstance(relative, str) or not relative:
        raise ProvenanceAuditError(f"missing {label} path")
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise ProvenanceAuditError(f"unsafe {label} path")
    try:
        path = (root / candidate).resolve(strict=True)
        path.relative_to(root)
    except (OSError, ValueError) as error:
        raise ProvenanceAuditError(f"unavailable {label} path") from error
    return path


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != STATUS
        or config.get("dataset_id") != "gse281364"
        or config.get("model_id") != MODEL_ID
        or config.get("task_id") != "signed_hepg2_reporter_activity_transfer"
        or config.get("outcome_access_authorized") is not False
        or config.get("reporter_count_access_authorized") is not False
        or config.get("metric_calculation_authorized") is not False
        or config.get("model_fit_authorized") is not False
        or config.get("prediction_generation_authorized") is not False
        or config.get("shortlist_authorized") is not False
        or config.get("complementarity_authorized") is not False
        or config.get("promotion_authorized") is not False
        or config.get("terminal_disposition")
        != "blocked_official_weights_only_scanner_unsupported_POP_opcode"
    ):
        raise ProvenanceAuditError("authorization firewall differs")
    rekey = config.get("rekeyed_prediction", {})
    if (
        rekey.get("selected_elements") != 1033
        or rekey.get("source_locus_groups") != 1033
        or rekey.get("long_range_blocks") != 239
        or rekey.get("independent_fitted_predictions") != 1
        or rekey.get("schema_seed_repeats_are_independent") is not False
    ):
        raise ProvenanceAuditError("re-keyed prediction contract differs")
    future = config.get("future_outcome_authority", {})
    if (
        future.get("role") != "exposed_development_MPRA_only"
        or future.get("read_during_this_audit") is not False
    ):
        raise ProvenanceAuditError("future outcome firewall differs")


def validate_archive_record(
    record: Mapping[str, Any], config: Mapping[str, Any]
) -> dict[str, Any]:
    authority = config["zenodo_record"]
    if record.get("id") != authority["record_id"]:
        raise ProvenanceAuditError("Zenodo record identifier differs")
    license_value = record.get("metadata", {}).get("license", {})
    if license_value.get("id") != authority["record_license"]:
        raise ProvenanceAuditError("Zenodo record license differs")
    files = {item.get("key"): item for item in record.get("files", ())}
    normalized: dict[str, Any] = {}
    for key in ("final_dump", "code_archive"):
        expected = authority[key]
        observed = files.get(expected["filename"], {})
        if (
            observed.get("size") != expected["size_bytes"]
            or observed.get("checksum") != f"md5:{expected['md5']}"
            or observed.get("links", {}).get("self") != expected["url"]
        ):
            raise ProvenanceAuditError(f"Zenodo {key} record differs")
        normalized[key] = {
            "filename": expected["filename"],
            "size_bytes": observed["size"],
            "md5": expected["md5"],
            "content_url": observed["links"]["self"],
        }
    return {
        "record_id": record["id"],
        "record_url": authority["record_url"],
        "record_created": record.get("created"),
        "record_updated": record.get("updated"),
        "record_license": authority["record_license"],
        "files": normalized,
    }


def validate_archive_file(path: Path, authority: Mapping[str, Any], *, label: str) -> None:
    if (
        not path.is_file()
        or path.is_symlink()
        or path.stat().st_size != authority["size_bytes"]
        or file_digest(path, "md5") != authority["md5"]
    ):
        raise ProvenanceAuditError(f"{label} bytes differ")


def validate_zip_member(
    archive: zipfile.ZipFile,
    *,
    name: str,
    size: int,
    compressed_size: int,
    crc32: str,
) -> zipfile.ZipInfo:
    try:
        member = archive.getinfo(name)
    except KeyError as error:
        raise ProvenanceAuditError(f"missing archive member: {name}") from error
    if (
        member.file_size != size
        or member.compress_size != compressed_size
        or f"{member.CRC:08x}" != crc32
        or member.is_dir()
    ):
        raise ProvenanceAuditError(f"archive member identity differs: {name}")
    return member


def read_archives(
    final_dump: Path,
    code_archive: Path,
    config: Mapping[str, Any],
) -> tuple[bytes, bytes, dict[str, bytes], list[dict[str, Any]]]:
    final_authority = config["zenodo_record"]["final_dump"]
    code_authority = config["zenodo_record"]["code_archive"]
    validate_archive_file(final_dump, final_authority, label="Zenodo final_dump")
    validate_archive_file(code_archive, code_authority, label="Zenodo code archive")
    inventory: list[dict[str, Any]] = []
    with zipfile.ZipFile(final_dump) as archive:
        if len(archive.infolist()) != final_authority["member_count"]:
            raise ProvenanceAuditError("final_dump member census differs")
        checkpoint_info = validate_zip_member(
            archive,
            name=final_authority["checkpoint_member"],
            size=final_authority["checkpoint_size_bytes"],
            compressed_size=final_authority["checkpoint_compressed_size_bytes"],
            crc32=final_authority["checkpoint_crc32"],
        )
        config_info = validate_zip_member(
            archive,
            name=final_authority["config_member"],
            size=final_authority["config_size_bytes"],
            compressed_size=final_authority["config_compressed_size_bytes"],
            crc32=final_authority["config_crc32"],
        )
        author_checkpoint = archive.read(checkpoint_info)
        author_config = archive.read(config_info)
        for member, role in (
            (checkpoint_info, "author_checkpoint"),
            (config_info, "author_training_config"),
        ):
            inventory.append(
                {
                    "archive": final_authority["filename"],
                    "member": member.filename,
                    "role": role,
                    "size_bytes": member.file_size,
                    "compressed_size_bytes": member.compress_size,
                    "crc32": f"{member.CRC:08x}",
                    "sha256": bytes_sha256(
                        author_checkpoint if role == "author_checkpoint" else author_config
                    ),
                }
            )
    code_values: dict[str, bytes] = {}
    with zipfile.ZipFile(code_archive) as archive:
        if len(archive.infolist()) != code_authority["member_count"]:
            raise ProvenanceAuditError("code-archive member census differs")
        for name, expected_hash in config["code_members"].items():
            try:
                member = archive.getinfo(name)
                value = archive.read(member)
            except KeyError as error:
                raise ProvenanceAuditError(f"missing code member: {name}") from error
            observed_hash = bytes_sha256(value)
            if observed_hash != expected_hash:
                raise ProvenanceAuditError(f"code member hash differs: {name}")
            code_values[name] = value
            inventory.append(
                {
                    "archive": code_authority["filename"],
                    "member": name,
                    "role": "author_code_or_environment",
                    "size_bytes": member.file_size,
                    "compressed_size_bytes": member.compress_size,
                    "crc32": f"{member.CRC:08x}",
                    "sha256": observed_hash,
                }
            )
    return author_checkpoint, author_config, code_values, inventory


def inspect_pickle_program(
    data: bytes,
    *,
    declared_globals: Sequence[str],
    allowed_opcodes: Sequence[str],
) -> dict[str, Any]:
    try:
        operations = list(pickletools.genops(data))
    except ValueError as error:
        raise ProvenanceAuditError("checkpoint pickle cannot be disassembled") from error
    observed_opcodes = {operation.name for operation, _, _ in operations}
    observed_globals = {
        str(argument)
        for operation, argument, _ in operations
        if operation.name == "GLOBAL"
    }
    if observed_opcodes != set(allowed_opcodes):
        raise ProvenanceAuditError("checkpoint pickle opcode set differs")
    if observed_globals != set(declared_globals):
        raise ProvenanceAuditError("checkpoint pickle global set differs")
    getattr_indices = [
        index
        for index, (operation, argument, _) in enumerate(operations)
        if operation.name == "GLOBAL" and argument == "__builtin__ getattr"
    ]
    if len(getattr_indices) != 1:
        raise ProvenanceAuditError("checkpoint getattr global census differs")
    index = getattr_indices[0]
    prior_strings = {
        str(argument)
        for operation, argument, _ in operations[max(0, index - 20) : index]
        if operation.name == "BINUNICODE"
    }
    next_globals = [
        str(argument)
        for operation, argument, _ in operations[index + 1 : index + 8]
        if operation.name == "GLOBAL"
    ]
    all_strings = {
        str(argument)
        for operation, argument, _ in operations
        if operation.name == "BINUNICODE"
    }
    if (
        "anneal_func" not in prior_strings
        or not next_globals
        or next_globals[0] != "torch.optim.lr_scheduler OneCycleLR"
        or "_annealing_cos" not in all_strings
    ):
        raise ProvenanceAuditError("builtins.getattr is not confined to OneCycleLR state")
    opcode_fingerprint = bytes_sha256(
        "\n".join(
            f"{operation.name}\t{argument!r}"
            for operation, argument, _ in operations
        ).encode("utf-8")
    )
    global_fingerprint = bytes_sha256(
        "\n".join(sorted(observed_globals)).encode("utf-8")
    )
    return {
        "pickle_protocol": int(data[1]) if data[:1] == b"\x80" else None,
        "opcode_count": len(operations),
        "opcode_counts": dict(
            sorted(Counter(operation.name for operation, _, _ in operations).items())
        ),
        "observed_opcodes": sorted(observed_opcodes),
        "declared_globals": sorted(observed_globals),
        "opcode_fingerprint_sha256": opcode_fingerprint,
        "global_fingerprint_sha256": global_fingerprint,
        "builtins_getattr_occurrences": 1,
        "builtins_getattr_scope": "OneCycleLR.anneal_func_to__annealing_cos",
    }


def inspect_torch_checkpoint_serialization(
    checkpoint: bytes, contract: Mapping[str, Any]
) -> dict[str, Any]:
    if bytes_sha256(checkpoint) != contract["checkpoint_sha256"]:
        raise ProvenanceAuditError("author checkpoint SHA-256 differs")
    try:
        with zipfile.ZipFile(io.BytesIO(checkpoint)) as archive:
            infos = archive.infolist()
            names = {member.filename for member in infos}
            validate_inner_member_names(
                names,
                contract["inner_tensor_storage_members"]
            )
            if (
                len(infos) != contract["inner_member_count"]
                or any(member.is_dir() or member.flag_bits & 0x1 for member in infos)
                or archive.testzip() is not None
            ):
                raise ProvenanceAuditError("author checkpoint inner archive differs")
            inventory = "\n".join(
                f"{member.filename}\t{member.file_size}\t{member.compress_size}\t{member.CRC:08x}"
                for member in infos
            ).encode("utf-8")
            data = archive.read("archive/data.pkl")
            version = archive.read("archive/version")
    except (OSError, zipfile.BadZipFile, RuntimeError) as error:
        raise ProvenanceAuditError("author checkpoint ZIP is invalid") from error
    if (
        bytes_sha256(inventory) != contract["inner_inventory_sha256"]
        or bytes_sha256(data) != contract["data_pickle_sha256"]
        or version.decode("ascii").strip() != contract["torch_archive_version"]
    ):
        raise ProvenanceAuditError("author checkpoint inner member bytes differ")
    result = inspect_pickle_program(
        data,
        declared_globals=contract["declared_globals"],
        allowed_opcodes=contract["allowed_opcodes"],
    )
    if (
        result["pickle_protocol"] != contract["pickle_protocol"]
        or result["opcode_fingerprint_sha256"]
        != contract["opcode_fingerprint_sha256"]
        or result["global_fingerprint_sha256"]
        != contract["global_fingerprint_sha256"]
    ):
        raise ProvenanceAuditError("author checkpoint pickle fingerprint differs")
    return {
        **result,
        "inner_member_count": len(infos),
        "tensor_storage_member_count": contract["inner_tensor_storage_members"],
        "inner_inventory_sha256": bytes_sha256(inventory),
        "data_pickle_sha256": bytes_sha256(data),
        "unexpected_archive_members": 0,
        "undeclared_globals": 0,
        "undeclared_opcodes": 0,
    }


def expected_inner_member_names(storage_members: int) -> set[str]:
    if storage_members < 1:
        raise ProvenanceAuditError("checkpoint storage-member census differs")
    return {
        "archive/data.pkl",
        "archive/version",
        *{f"archive/data/{index}" for index in range(storage_members)},
    }


def validate_inner_member_names(
    observed_names: set[str], storage_members: int
) -> None:
    if observed_names != expected_inner_member_names(storage_members):
        raise ProvenanceAuditError("author checkpoint inner member names differ")


def validate_architecture_config(
    author_config_bytes: bytes, hf_config: Mapping[str, Any]
) -> dict[str, Any]:
    try:
        author = json.loads(author_config_bytes)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ProvenanceAuditError("invalid author training config") from error
    architecture = (
        "stem_ch",
        "stem_ks",
        "ef_ks",
        "ef_block_sizes",
        "resize_factor",
        "pool_sizes",
    )
    if (
        any(author.get(key) != hf_config.get(key) for key in architecture)
        or author.get("use_reverse_channel") is not False
        or author.get("reverse_augment") is not True
        or author.get("use_shift") is not True
        or hf_config.get("in_ch") != 4
        or hf_config.get("activation") != "SiLU"
    ):
        raise ProvenanceAuditError("author and repackaged architecture config differ")
    return {
        "architecture_fields_exact": True,
        "author_reverse_augment": True,
        "author_training_shift_augment": True,
        "author_reverse_channel": False,
        "input_channels": 4,
        "activation": "SiLU",
    }


@contextmanager
def network_guard():
    original_socket = socket.socket
    original_connection = socket.create_connection
    attempts: list[str] = []

    class GuardedSocket(original_socket):
        def __new__(cls, family=-1, *args, **kwargs):
            if family in {socket.AF_INET, socket.AF_INET6}:
                attempts.append(str(family))
                raise ProvenanceAuditError("network access during checkpoint conversion")
            return super().__new__(cls, family, *args, **kwargs)

    def blocked_connection(*args, **kwargs):
        del args, kwargs
        attempts.append("create_connection")
        raise ProvenanceAuditError("network access during checkpoint conversion")

    socket.socket = GuardedSocket
    socket.create_connection = blocked_connection
    try:
        yield attempts
    finally:
        socket.socket = original_socket
        socket.create_connection = original_connection


def _safe_global_objects() -> dict[str, object]:
    return {
        "builtins.getattr": builtins.getattr,
        "collections.defaultdict": defaultdict,
        "torch.optim.adamw.AdamW": torch.optim.AdamW,
        "torch.optim.lr_scheduler.OneCycleLR": torch.optim.lr_scheduler.OneCycleLR,
    }


def probe_official_unsafe_global_scanner(checkpoint: bytes) -> dict[str, Any]:
    try:
        observed = sorted(
            torch.serialization.get_unsafe_globals_in_checkpoint(
                io.BytesIO(checkpoint)
            )
        )
    except pickle.UnpicklingError as error:
        if str(error) != "Unsupported operand 48":
            raise ProvenanceAuditError(
                "official weights-only scanner failed unexpectedly"
            ) from error
        return {
            "status": "blocked_official_weights_only_scanner_unsupported_POP_opcode",
            "official_scanner": (
                "torch.serialization.get_unsafe_globals_in_checkpoint"
            ),
            "exception_type": "_pickle.UnpicklingError",
            "exception_message": str(error),
            "unsupported_opcode": "POP",
            "unsupported_opcode_decimal": 48,
            "torch_reported_unsafe_globals": None,
            "torch_load_invoked": False,
            "weights_only_load_attempted": False,
            "weights_only_false_used": False,
            "custom_deserializer_used": False,
            "unrestricted_deserializer_used": False,
        }
    return {
        "status": "pass_official_weights_only_scanner",
        "official_scanner": "torch.serialization.get_unsafe_globals_in_checkpoint",
        "torch_reported_unsafe_globals": observed,
        "torch_load_invoked": False,
        "weights_only_load_attempted": False,
        "weights_only_false_used": False,
        "custom_deserializer_used": False,
        "unrestricted_deserializer_used": False,
    }


def _extract_tensor_state(value: object) -> dict[str, torch.Tensor]:
    if not isinstance(value, dict) or not isinstance(value.get("state_dict"), dict):
        raise ProvenanceAuditError("author checkpoint state_dict differs")
    state = value["state_dict"]
    if any(
        not isinstance(key, str) or not isinstance(tensor, torch.Tensor)
        for key, tensor in state.items()
    ):
        raise ProvenanceAuditError("author checkpoint state_dict has non-tensor leaves")
    return dict(state)


def _state_dict_from_author_checkpoint(
    checkpoint: bytes,
    *,
    expected_unsafe_globals: Sequence[str] = (),
    safe_global_names: Sequence[str] = (),
) -> tuple[dict[str, torch.Tensor], dict[str, Any]]:
    observed_unsafe = sorted(
        torch.serialization.get_unsafe_globals_in_checkpoint(io.BytesIO(checkpoint))
    )
    if observed_unsafe != sorted(expected_unsafe_globals):
        raise ProvenanceAuditError("torch unsafe-global inventory differs")
    available = _safe_global_objects()
    if set(safe_global_names) != set(expected_unsafe_globals) or any(
        name not in available for name in safe_global_names
    ):
        raise ProvenanceAuditError("narrow safe-global allowlist differs")
    objects = [available[name] for name in safe_global_names]
    try:
        with network_guard() as first_attempts:
            with torch.serialization.safe_globals(objects):
                first = torch.load(
                    io.BytesIO(checkpoint), map_location="cpu", weights_only=True
                )
        with network_guard() as second_attempts:
            with torch.serialization.safe_globals(objects):
                second = torch.load(
                    io.BytesIO(checkpoint), map_location="cpu", weights_only=True
                )
    except Exception as error:
        raise ProvenanceAuditError(
            "author checkpoint failed narrow allowlisted weights-only loading"
        ) from error
    first_state = _extract_tensor_state(first)
    second_state = _extract_tensor_state(second)
    deterministic, _, summary = compare_checkpoint_states(first_state, second_state)
    if not deterministic or first_attempts or second_attempts:
        raise ProvenanceAuditError("checkpoint conversion is not deterministic or offline")
    return first_state, {
        "torch_reported_unsafe_globals": observed_unsafe,
        "narrow_safe_globals": sorted(safe_global_names),
        "weights_only": True,
        "weights_only_false_used": False,
        "network_attempts": 0,
        "repeat_conversion_bitwise_equal": True,
        "repeat_canonical_state_sha256": summary[
            "author_canonical_state_sha256"
        ],
    }


def tensor_digest(tensor: torch.Tensor) -> str:
    return bytes_sha256(save({"tensor": tensor.detach().cpu().contiguous()}))


def compare_checkpoint_states(
    author_state: Mapping[str, torch.Tensor], hf_state: Mapping[str, torch.Tensor]
) -> tuple[bool, list[dict[str, Any]], dict[str, Any]]:
    names = sorted(set(author_state).union(hf_state))
    rows: list[dict[str, Any]] = []
    all_equal = set(author_state) == set(hf_state)
    total_elements = 0
    for name in names:
        author = author_state.get(name)
        repackaged = hf_state.get(name)
        present = author is not None and repackaged is not None
        shape_equal = present and tuple(author.shape) == tuple(repackaged.shape)
        dtype_equal = present and author.dtype == repackaged.dtype
        value_equal = bool(
            present and shape_equal and dtype_equal and torch.equal(author, repackaged)
        )
        all_equal = all_equal and value_equal
        if author is not None:
            total_elements += author.numel()
        rows.append(
            {
                "tensor_name": name,
                "author_present": str(author is not None).lower(),
                "repackaged_present": str(repackaged is not None).lower(),
                "shape": "x".join(map(str, author.shape)) if author is not None else "missing",
                "dtype": str(author.dtype) if author is not None else "missing",
                "author_tensor_sha256": tensor_digest(author) if author is not None else "missing",
                "repackaged_tensor_sha256": (
                    tensor_digest(repackaged) if repackaged is not None else "missing"
                ),
                "bitwise_equal": str(value_equal).lower(),
            }
        )
    ordered_author = {
        key: author_state[key].detach().cpu().contiguous() for key in sorted(author_state)
    }
    ordered_hf = {key: hf_state[key].detach().cpu().contiguous() for key in sorted(hf_state)}
    summary = {
        "tensor_names_exact": set(author_state) == set(hf_state),
        "tensor_count": len(author_state),
        "parameter_and_buffer_elements": total_elements,
        "all_tensors_bitwise_equal": all_equal,
        "author_canonical_state_sha256": bytes_sha256(save(ordered_author)),
        "repackaged_canonical_state_sha256": bytes_sha256(save(ordered_hf)),
        "restricted_weights_only_load": True,
        "unsafe_pickle_fallback_used": False,
    }
    return all_equal, rows, summary


def validate_project_authorities(
    project_root: Path, config: Mapping[str, Any]
) -> dict[str, Path]:
    trees: dict[str, Path] = {}
    for label, authority in config["immutable_trees"].items():
        tree = safe_project_path(project_root, authority["path"], label=label)
        if file_digest(tree / "ARTIFACTS.json") != authority["artifacts_sha256"]:
            raise ProvenanceAuditError(f"immutable tree manifest differs: {label}")
        verify_frozen_tree(tree)
        trees[label] = tree
    for relative, expected in config["execution_sources"].items():
        path = safe_project_path(project_root, relative, label="execution source")
        if file_digest(path) != expected:
            raise ProvenanceAuditError(f"execution source differs: {relative}")
    admission = load_json(
        trees["admission"] / "validation/admission_receipt.json",
        label="admission receipt",
    )
    fixture = load_json(trees["fixture"] / "fixture/receipt.json", label="fixture receipt")
    prediction = load_json(
        trees["predictions"] / "predictions/receipt.json", label="prediction receipt"
    )
    reconciliation = load_json(
        trees["reconciliation"] / "reconciliation/receipt.json",
        label="reconciliation receipt",
    )
    rekey = config["rekeyed_prediction"]
    selected = trees["reconciliation"] / rekey["member"]
    if (
        file_digest(selected) != rekey["member_sha256"]
        or admission.get("checkpoint_sha256") != HF_CHECKPOINT_SHA256
        or not admission.get("author_and_hf_model_source_byte_identical")
        or admission.get("observed_outcomes_loaded")
        or admission.get("sealed_outcomes_loaded")
        or fixture.get("elements") != 4359
        or fixture.get("outcomes_read")
        or fixture.get("reporter_counts_read")
        or prediction.get("checkpoint_sha256") != HF_CHECKPOINT_SHA256
        or prediction.get("elements") != 4359
        or prediction.get("outcomes_read")
        or prediction.get("reporter_counts_read")
        or prediction.get("sealed_outcomes_read")
        or prediction.get("model_fitted_or_adapted")
        or reconciliation.get("selected_elements") != 1033
        or reconciliation.get("row_keys") != 10330
        or reconciliation.get("outcomes_read")
        or reconciliation.get("metrics_calculated")
        or reconciliation.get("model_fit")
        or reconciliation.get("new_predictions_generated")
    ):
        raise ProvenanceAuditError("project execution or outcome firewall differs")
    return trees


def validate_preserved_failure_reuse(
    project_root: Path,
    config: Mapping[str, Any],
    *,
    record: Path,
    final_dump: Path,
    code_archive: Path,
) -> dict[str, Any]:
    contract = config.get("preserved_failure_reuse", {})
    if (
        contract.get("authorized") is not True
        or contract.get("source_job_id") != 21091805
        or contract.get("source_state") != "FAILED"
        or contract.get("source_exit_code") != "1:0"
        or contract.get("failure_reason")
        != "restricted_weights_only_loader_rejected_undeclared_builtins_getattr"
        or contract.get("reuse_requires_reverification") is not True
        or contract.get("reuse_is_not_a_frozen_result_authority") is not True
    ):
        raise ProvenanceAuditError("preserved-failure reuse is not authorized")
    source = safe_project_path(
        project_root, contract["source_path"], label="preserved failure"
    )
    expected = {
        "record": source / contract["record_member"],
        "final_dump": source / contract["final_dump_member"],
        "code_archive": source / contract["code_archive_member"],
    }
    observed = {
        "record": record.resolve(strict=True),
        "final_dump": final_dump.resolve(strict=True),
        "code_archive": code_archive.resolve(strict=True),
    }
    if observed != expected:
        raise ProvenanceAuditError("retry inputs do not resolve to preserved failure")
    stderr = source / contract["audit_stderr_member"]
    if (
        file_digest(stderr) != contract["audit_stderr_sha256"]
        or file_digest(observed["record"])
        != contract["record_sha256"]
        or observed["final_dump"].stat().st_size
        != contract["final_dump_size_bytes"]
        or file_digest(observed["final_dump"], "md5")
        != contract["final_dump_md5"]
        or observed["code_archive"].stat().st_size
        != contract["code_archive_size_bytes"]
        or file_digest(observed["code_archive"], "md5")
        != contract["code_archive_md5"]
    ):
        raise ProvenanceAuditError("preserved-failure bytes differ")
    return {
        "schema_version": "masld-bench-preserved-failure-reuse-contract-v1",
        "status": "pass_exact_failed_download_reuse",
        "source_path": contract["source_path"],
        "source_job_id": contract["source_job_id"],
        "source_state": contract["source_state"],
        "source_exit_code": contract["source_exit_code"],
        "failure_reason": contract["failure_reason"],
        "record_sha256": contract["record_sha256"],
        "final_dump_size_bytes": contract["final_dump_size_bytes"],
        "final_dump_md5": contract["final_dump_md5"],
        "code_archive_size_bytes": contract["code_archive_size_bytes"],
        "code_archive_md5": contract["code_archive_md5"],
        "redownloaded": False,
        "reverified_before_reuse": True,
        "source_is_result_authority": False,
    }


def validate_terminal_scanner_failure(
    project_root: Path, config: Mapping[str, Any]
) -> dict[str, Any]:
    contract = config.get("terminal_scanner_failure", {})
    if (
        contract.get("source_job_id") != 21092032
        or contract.get("source_state") != "FAILED"
        or contract.get("source_exit_code") != "1:0"
        or contract.get("official_scanner")
        != "torch.serialization.get_unsafe_globals_in_checkpoint"
        or contract.get("exception_type") != "_pickle.UnpicklingError"
        or contract.get("exception_message") != "Unsupported operand 48"
        or contract.get("unsupported_opcode") != "POP"
        or contract.get("unsupported_opcode_decimal") != 48
        or contract.get("torch_load_reached") is not False
        or contract.get("custom_deserializer_allowed") is not False
        or contract.get("unrestricted_deserializer_allowed") is not False
    ):
        raise ProvenanceAuditError("terminal scanner-failure contract differs")
    source = safe_project_path(
        project_root, contract["source_path"], label="scanner failure"
    )
    stderr = source / contract["audit_stderr_member"]
    tests = source / contract["tests_stderr_member"]
    if (
        file_digest(stderr) != contract["audit_stderr_sha256"]
        or file_digest(tests) != contract["tests_stderr_sha256"]
        or "Unsupported operand 48" not in stderr.read_text(encoding="utf-8")
        or "Ran 12 tests" not in tests.read_text(encoding="utf-8")
        or "OK" not in tests.read_text(encoding="utf-8")
    ):
        raise ProvenanceAuditError("terminal scanner-failure evidence differs")
    return {
        "schema_version": "masld-bench-terminal-scanner-failure-evidence-v1",
        "status": "pass_exact_preserved_scanner_failure_evidence",
        **dict(contract),
        "source_is_result_authority": False,
    }


def build_task_spec(
    config: Mapping[str, Any],
    *,
    author_checkpoint_sha256: str,
    canonical_state_sha256: str,
) -> dict[str, Any]:
    contract = config["pre_scoring_contract"]
    return {
        "schema_version": "masld-bench-pre-scoring-task-spec-v1",
        "status": "frozen_descriptive_pre_scoring_contract",
        "dataset_id": "gse281364",
        "model_id": MODEL_ID,
        "task_id": "signed_hepg2_reporter_activity_transfer",
        "scientific_role": "exposed_development_descriptive_native_output_transfer",
        "prediction_universe": {
            "elements": 1033,
            "source_locus_groups": 1033,
            "long_range_blocks": 239,
            "prediction_member": config["rekeyed_prediction"]["member"],
            "prediction_member_sha256": config["rekeyed_prediction"]["member_sha256"],
            "prediction_is_static_across_contexts": True,
            "independent_fitted_predictions": 1,
            "schema_seed_repeats_are_independent": False,
        },
        "checkpoint": {
            "author_record": "Zenodo:10558183",
            "author_member": AUTHOR_MEMBER,
            "author_checkpoint_sha256": author_checkpoint_sha256,
            "repackaged_checkpoint_sha256": HF_CHECKPOINT_SHA256,
            "canonical_state_sha256": canonical_state_sha256,
            "equivalence": "all_state_dict_tensors_bitwise_equal",
        },
        "endpoint": {
            "prediction": contract["prediction_endpoint"],
            "observed": contract["observed_endpoint"],
            "allele_sign": "ALT_minus_REF",
            "eligible_contexts": contract["eligible_contexts"],
            "excluded_contexts": contract["excluded_contexts"],
            "pseudocount": contract["pseudocount"],
        },
        "inference": {
            "unit": "element",
            "resampling_unit": contract["resampling_unit"],
            "biological_donors": contract["biological_donors"],
            "experimental_replicates_per_context": contract[
                "experimental_replicates_per_context"
            ],
            "replicates_are_independent_donors": contract[
                "replicates_are_independent_donors"
            ],
            "bootstrap_replicates": contract["bootstrap_replicates"],
            "bootstrap_seed": contract["bootstrap_seed"],
        },
        "descriptors": {
            "primary": contract["primary_descriptor"],
            "secondary": contract["secondary_descriptor"],
            "context_specific_only": True,
            "p_values_authorized": False,
            "multiple_testing": "not_applicable_no_hypothesis_tests_or_p_values",
        },
        "claims": {
            "models_ranked": False,
            "shortlist_eligible": False,
            "complementarity_eligible": False,
            "champion_eligible": False,
            "standalone_signed_eqtl_or_ieqtl_eligible": False,
            "treatment_interaction_prediction": False,
        },
    }


def build_evaluator_contract(config: Mapping[str, Any]) -> dict[str, Any]:
    future = config["future_outcome_authority"]
    return {
        "schema_version": "masld-bench-descriptive-evaluator-contract-v1",
        "status": "frozen_not_executed",
        "task_id": "signed_hepg2_reporter_activity_transfer",
        "prediction_input": {
            "tree": config["immutable_trees"]["reconciliation"]["path"],
            "tree_artifacts_sha256": config["immutable_trees"]["reconciliation"][
                "artifacts_sha256"
            ],
            "member": config["rekeyed_prediction"]["member"],
            "member_sha256": config["rekeyed_prediction"]["member_sha256"],
            "score_field": "alt_minus_ref_score",
            "deduplicate_to": ["element_id", "source_locus_group_id"],
            "expected_rows": 1033,
        },
        "outcome_input": {
            "tree": future["path"],
            "tree_artifacts_sha256": future["artifacts_sha256"],
            "role": future["role"],
            "read_during_provenance_audit": False,
        },
        "outcome_transform": {
            "replicate_delta": (
                "log2((ALT_RNA+0.5)/(ALT_DNA+0.5))-"
                "log2((REF_RNA+0.5)/(REF_DNA+0.5))"
            ),
            "aggregate": "arithmetic_mean_across_four_complete_assay_replicates",
            "missingness": "require_complete_REF_ALT_pairs_never_encode_missing_as_zero",
        },
        "evaluation": {
            "contexts": ["HepG2_control", "HepG2_PAOA"],
            "primary_descriptor": "spearman_rank_correlation",
            "secondary_descriptor": "pearson_correlation",
            "uncertainty": "10000_source_locus_group_bootstrap_replicates",
            "bootstrap_seed": 20260825,
            "outer_fold_results": "diagnostic_only_not_independent_replications",
            "schema_seed_rows": "deduplicated_not_independent_replications",
            "model_ranking": False,
            "promotion_gate": "none_descriptive_nonchampion",
        },
        "execution": {
            "authorized_in_this_audit": False,
            "outcomes_read": False,
            "metrics_calculated": False,
        },
    }


def audit(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise ProvenanceAuditError("audit output exists")
    config = load_config(arguments.config)
    validate_config(config)
    project_root = arguments.project_root.resolve(strict=True)
    trees = validate_project_authorities(project_root, config)
    reuse_contract = validate_preserved_failure_reuse(
        project_root,
        config,
        record=arguments.zenodo_record,
        final_dump=arguments.final_dump,
        code_archive=arguments.code_archive,
    )
    scanner_failure_evidence = validate_terminal_scanner_failure(
        project_root, config
    )
    record = load_json(arguments.zenodo_record, label="Zenodo record")
    normalized_record = validate_archive_record(record, config)
    author_checkpoint, author_config, code_values, archive_inventory = read_archives(
        arguments.final_dump, arguments.code_archive, config
    )
    serialization_contract = config["zenodo_record"]["final_dump"][
        "restricted_serialization"
    ]
    serialization_contract = {
        **serialization_contract,
        "checkpoint_sha256": config["zenodo_record"]["final_dump"][
            "checkpoint_sha256"
        ],
    }
    serialization_inventory = inspect_torch_checkpoint_serialization(
        author_checkpoint, serialization_contract
    )
    admission_model = (trees["admission"] / "sources/author_model.py").read_bytes()
    if code_values["human_legnet-main/model.py"] != admission_model:
        raise ProvenanceAuditError("executed model source differs from Zenodo author code")
    hf_config = load_json(
        trees["admission"] / "metadata/config.json", label="repackaged config"
    )
    config_equivalence = validate_architecture_config(author_config, hf_config)
    scanner_probe = probe_official_unsafe_global_scanner(author_checkpoint)
    if scanner_probe["status"] != config["terminal_disposition"]:
        raise ProvenanceAuditError(
            "official scanner no longer reproduces the terminal disposition"
        )
    provenance_pass = False
    tensors_equal = None
    tensor_summary = {
        "comparison_performed": False,
        "reason": config["terminal_disposition"],
        "author_tensor_state_loaded": False,
        "repackaged_tensor_state_loaded": False,
    }

    arguments.output.mkdir(mode=0o750)
    provenance = arguments.output / "provenance"
    provenance.mkdir(mode=0o750)
    write_json(provenance / "preserved_failure_reuse_contract.json", reuse_contract)
    write_json(
        provenance / "terminal_scanner_failure_evidence.json",
        scanner_failure_evidence,
    )
    write_json(provenance / "zenodo_record.json", normalized_record)
    write_json(
        provenance / "restricted_serialization_inventory.json",
        {
            **serialization_inventory,
            **scanner_probe,
            "narrow_safe_globals_considered_not_used": sorted(
                serialization_contract["narrow_safe_global_objects"]
            ),
        },
    )
    (provenance / "author_best_model_test1_val2.ckpt").write_bytes(author_checkpoint)
    (provenance / "author_training_config.json").write_bytes(author_config)
    write_tsv(
        provenance / "archive_member_inventory.tsv",
        tuple(archive_inventory[0]),
        archive_inventory,
    )
    author_checkpoint_hash = bytes_sha256(author_checkpoint)
    receipt = {
        "schema_version": SCHEMA,
        "status": config["terminal_disposition"],
        "dataset_id": "gse281364",
        "model_id": MODEL_ID,
        "zenodo_record_id": 10558183,
        "zenodo_final_dump_md5_verified": True,
        "zenodo_code_archive_md5_verified": True,
        "author_checkpoint_member": AUTHOR_MEMBER,
        "author_checkpoint_sha256": author_checkpoint_hash,
        "repackaged_checkpoint_sha256": HF_CHECKPOINT_SHA256,
        "checkpoint_packaging_byte_identical": False,
        "checkpoint_state_tensors_bitwise_equal": tensors_equal,
        "author_and_executed_model_source_byte_identical": True,
        "architecture_config_equivalent": config_equivalence[
            "architecture_fields_exact"
        ],
        "input_preprocessing_source_hashes_verified": True,
        "project_execution_source_hashes_verified": True,
        "restricted_serialization_inventory_passed": True,
        "torch_reported_unsafe_globals": None,
        "narrow_safe_globals_considered_not_used": sorted(
            serialization_contract["narrow_safe_global_objects"]
        ),
        "builtins_getattr_reason": serialization_contract["getattr_reason"],
        "official_scanner": scanner_probe["official_scanner"],
        "official_scanner_exception_type": scanner_probe["exception_type"],
        "official_scanner_exception_message": scanner_probe["exception_message"],
        "official_scanner_unsupported_opcode": scanner_probe[
            "unsupported_opcode"
        ],
        "torch_load_invoked": False,
        "weights_only_load_attempted": False,
        "weights_only_false_used": False,
        "custom_deserializer_used": False,
        "unrestricted_deserializer_used": False,
        "network_attempts_during_conversion": 0,
        "repeat_conversion_bitwise_equal": None,
        "training_domain": "HepG2_lentiMPRA",
        "evaluation_input_lane": "GRCh38p14_230bp_genomic_window_transfer",
        "native_output": "uncalibrated_lentiMPRA_reporter_expression_score",
        "single_checkpoint_not_author_ensemble": True,
        "zero_shift_inference": True,
        "reverse_complement_averaged": True,
        "tensors": tensor_summary,
        "provenance_equivalence_pass": provenance_pass,
        "descriptive_native_output_scoring_eligible": False,
        "native_assay_input_claim": False,
        "models_ranked": False,
        "shortlist_eligible": False,
        "complementarity_eligible": False,
        "standalone_champion_eligible": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_assets_read": False,
        "metrics_calculated": False,
        "model_fit": False,
        "predictions_generated": False,
        "schema_seed_repeats_are_independent": False,
        "independent_fitted_predictions": 1,
        "future_outcome_authority_declared_not_read": True,
    }
    write_json(provenance / "receipt.json", receipt)
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--zenodo-record", type=Path, required=True)
    parser.add_argument("--final-dump", type=Path, required=True)
    parser.add_argument("--code-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    audit(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
