"""Hash-first, offline staging for Geneformer checkpoint bundles.

This module never downloads output files.  It accepts a user-populated source
tree, verifies every required object against the frozen upstream manifest, and
publishes only verified files.  Auxiliary pickle dictionaries are decoded with
all global/class construction disabled, validated as plain JSON-compatible
mappings, and re-serialized as both canonical JSON and project-owned pickle.
"""

from __future__ import annotations

import io
import json
import math
import os
from pathlib import Path
import pickle
import shutil
import struct
import tempfile
from typing import Any, Mapping

from .artifacts import (
    ArtifactError,
    freeze_tree,
    publish_directory_noreplace,
    reject_symlink_components,
    verify_frozen_tree,
    write_bytes_exclusive,
    write_json_exclusive,
)
from .hashing import canonical_sha256, require_sha256, sha256_bytes, sha256_file


class CheckpointPreflightError(RuntimeError):
    """Raised when a checkpoint source cannot be staged safely."""


_GENEFORMER_VARIANTS = {
    "geneformer_v1_10m": (0, "V1", None),
    "geneformer_v2_104m": (1, "V2_shared", "V2_104M"),
    "geneformer_v2_316m": (2, "V2_shared", "V2_316M"),
}


# Upstream Geneformer stores its gene-median and token dictionaries as plain
# str -> numeric mappings whose values were pickled as numpy scalars.  Decoding
# them needs exactly two globals and nothing else.  They are included under a
# closed allowlist with a hard rejection of object and void dtypes, which is
# the dtype family that makes numpy scalar reconstruction dangerous: an object
# dtype causes the buffer itself to be unpickled, reintroducing arbitrary
# construction.  Every included value is immediately coerced to a plain Python
# int or float, so no numpy object survives into a staged bundle.
# numpy renamed `numpy.core` to `numpy._core`, so the same reconstructor
# appears under two module paths depending on the numpy that wrote the pickle.
# Upstream Geneformer's dictionaries carry the legacy spelling; anything
# re-pickled under numpy 2.x carries the underscored one.  Both are the same
# function and both are subject to the identical dtype rejection below.
_ALLOWED_PICKLE_GLOBALS = frozenset(
    {
        ("numpy.core.multiarray", "scalar"),
        ("numpy._core.multiarray", "scalar"),
        ("numpy", "dtype"),
    }
)
# Numeric scalar reconstruction is emulated with `struct` rather than by
# importing numpy.  That keeps numpy out of the security-critical path
# entirely: no numpy code executes while decoding an untrusted pickle, and the
# stager keeps working under the bare project interpreter, which has no numpy.
# Only fixed-width numeric dtypes are representable here; object, void, string,
# and datetime dtypes have no entry and are therefore rejected by construction.
_NUMERIC_DTYPE_FORMATS = {
    "b1": "?",
    "f4": "f",
    "f8": "d",
    "i1": "b",
    "i2": "h",
    "i4": "i",
    "i8": "q",
    "u1": "B",
    "u2": "H",
    "u4": "I",
    "u8": "Q",
}
_BYTE_ORDERS = {"<": "<", ">": ">", "=": "=", "|": "="}


class _RestrictedDtype:
    """A numeric-only stand-in for numpy.dtype during restricted unpickling."""

    __slots__ = ("code", "byte_order")

    def __init__(self, spec: Any) -> None:
        if not isinstance(spec, str) or spec not in _NUMERIC_DTYPE_FORMATS:
            raise CheckpointPreflightError(
                f"pickle numpy dtype is not a permitted numeric type: {spec!r}"
            )
        self.code = spec
        self.byte_order = "="

    def __setstate__(self, state: Any) -> None:
        # numpy's dtype state is a tuple whose second element is byte order.
        if not isinstance(state, tuple) or len(state) < 2:
            raise CheckpointPreflightError("pickle numpy dtype state is invalid")
        order = state[1]
        if order not in _BYTE_ORDERS:
            raise CheckpointPreflightError(
                f"pickle numpy dtype byte order is invalid: {order!r}"
            )
        self.byte_order = _BYTE_ORDERS[order]


def _restricted_numpy_dtype(spec: Any = None, *args: Any, **kwargs: Any) -> Any:
    return _RestrictedDtype(spec)


def _restricted_numpy_scalar(dtype: Any, payload: Any) -> Any:
    if not isinstance(dtype, _RestrictedDtype):
        raise CheckpointPreflightError(
            "pickle numpy scalar requires a permitted numeric dtype"
        )
    if not isinstance(payload, (bytes, bytearray)):
        raise CheckpointPreflightError(
            "pickle numpy scalar payload must be raw bytes"
        )
    fmt = f"{dtype.byte_order}{_NUMERIC_DTYPE_FORMATS[dtype.code]}"
    expected = struct.calcsize(fmt)
    if len(payload) != expected:
        raise CheckpointPreflightError(
            f"pickle numpy scalar payload must be {expected} bytes"
        )
    (value,) = struct.unpack(fmt, bytes(payload))
    if isinstance(value, float) and not math.isfinite(value):
        raise CheckpointPreflightError("pickle numpy float is not finite")
    return value


class _PlainDataUnpickler(pickle.Unpickler):
    """Deny every pickle global except two numeric-scalar reconstructors."""

    def find_class(self, module: str, name: str) -> Any:
        if (module, name) not in _ALLOWED_PICKLE_GLOBALS:
            raise CheckpointPreflightError(
                f"pickle global/class construction is forbidden: {module}.{name}"
            )
        if name == "dtype":
            return _restricted_numpy_dtype
        return _restricted_numpy_scalar


def _plain_json_value(value: Any, *, path: str = "root") -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise CheckpointPreflightError(f"non-finite pickle value at {path}")
        return value
    if isinstance(value, dict):
        result: dict[str, Any] = {}
        for key, child in value.items():
            # Keys must be strings, but the EMPTY string is included verbatim.
            # Upstream Geneformer's V1 gene_name_id_dict_gc30M contains exactly
            # one such entry, '' -> ENSG00000285325, a real Ensembl gene whose
            # symbol is blank.  It is inert for tokenization, since no real gene
            # symbol can look it up, and JSON represents "" as a key without
            # difficulty.  Dropping it would silently alter a hash-verified
            # upstream output file, which is worse than carrying it.
            if not isinstance(key, str):
                raise CheckpointPreflightError(
                    f"pickle mapping keys must be strings at {path}"
                )
            result[key] = _plain_json_value(child, path=f"{path}.{key}")
        return result
    raise CheckpointPreflightError(
        f"pickle contains forbidden {type(value).__name__} at {path}"
    )


def load_restricted_plain_mapping(payload: bytes) -> dict[str, Any]:
    """Decode one plain dictionary while denying pickle globals and classes."""

    stream = io.BytesIO(payload)
    try:
        value = _PlainDataUnpickler(stream).load()
    except CheckpointPreflightError:
        raise
    except Exception as error:
        raise CheckpointPreflightError(
            f"auxiliary pickle cannot be decoded as restricted plain data: {error}"
        ) from error
    if stream.read(1):
        raise CheckpointPreflightError(
            "auxiliary pickle contains trailing content after the first object"
        )
    normalized = _plain_json_value(value)
    if not isinstance(normalized, dict) or not normalized:
        raise CheckpointPreflightError(
            "auxiliary pickle must contain one nonempty plain dictionary"
        )
    return normalized


def _load_manifest(path: Path) -> Mapping[str, Any]:
    source = reject_symlink_components(path, label="checkpoint manifest")
    try:
        raw = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CheckpointPreflightError(
            f"checkpoint manifest is unreadable: {error}"
        ) from error
    if not isinstance(raw, Mapping) or raw.get("schema_version") != (
        "masld-bench-upstream-checkpoint-preflight-v1"
    ):
        raise CheckpointPreflightError("unsupported Geneformer checkpoint manifest")
    if raw.get("repository") != "ctheodoris/Geneformer":
        raise CheckpointPreflightError("checkpoint manifest repository differs")
    revision = raw.get("repository_revision")
    if (
        not isinstance(revision, str)
        or len(revision) != 40
        or any(character not in "0123456789abcdef" for character in revision)
    ):
        raise CheckpointPreflightError("checkpoint repository revision is not exact")
    return raw


def _artifact_spec(raw: Any, *, label: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping) or {
        "path",
        "sha256",
        "size_bytes",
    }.difference(raw):
        raise CheckpointPreflightError(f"{label} artifact contract is incomplete")
    relative = Path(str(raw["path"]))
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise CheckpointPreflightError(f"{label} artifact path is unsafe")
    size = raw["size_bytes"]
    if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        raise CheckpointPreflightError(f"{label} artifact size is invalid")
    return {
        "path": relative.as_posix(),
        "sha256": require_sha256(raw["sha256"], field_name=f"{label}.sha256"),
        "size_bytes": size,
    }


def _required_artifacts(
    manifest: Mapping[str, Any], model_id: str
) -> tuple[dict[str, Any], ...]:
    try:
        weight_index, shared_key, variant_key = _GENEFORMER_VARIANTS[model_id]
    except KeyError as error:
        raise CheckpointPreflightError(
            f"unsupported Geneformer model_id: {model_id}"
        ) from error
    weights = manifest.get("artifacts")
    auxiliary = manifest.get("auxiliary_artifacts")
    if not isinstance(weights, list) or not isinstance(auxiliary, Mapping):
        raise CheckpointPreflightError("checkpoint manifest artifact inventory is invalid")
    try:
        selected: list[Any] = [weights[weight_index], *auxiliary[shared_key]]
        if variant_key is not None:
            selected.extend(auxiliary[variant_key])
    except (IndexError, KeyError, TypeError) as error:
        raise CheckpointPreflightError(
            f"checkpoint manifest lacks the {model_id} artifact roster"
        ) from error
    normalized = tuple(
        _artifact_spec(raw, label=f"{model_id}[{index}]")
        for index, raw in enumerate(selected)
    )
    paths = [item["path"] for item in normalized]
    if len(set(paths)) != len(paths):
        raise CheckpointPreflightError("checkpoint manifest repeats an artifact path")
    return normalized


def _verified_source(path: Path, spec: Mapping[str, Any]) -> Path:
    source = reject_symlink_components(path, label="checkpoint source artifact")
    if source.is_symlink() or not source.is_file():
        raise CheckpointPreflightError(f"checkpoint source is missing: {spec['path']}")
    before = source.stat(follow_symlinks=False)
    if before.st_size != spec["size_bytes"]:
        raise CheckpointPreflightError(
            f"checkpoint source size differs: {spec['path']}"
        )
    digest = sha256_file(source)
    after = source.stat(follow_symlinks=False)
    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if identity_before != identity_after or source.is_symlink():
        raise CheckpointPreflightError(
            f"checkpoint source changed while hashing: {spec['path']}"
        )
    if digest != spec["sha256"]:
        raise CheckpointPreflightError(
            f"checkpoint source SHA-256 differs: {spec['path']}"
        )
    return source


def _copy_exclusive(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    reject_symlink_components(target, label="checkpoint staged artifact")
    with source.open("rb") as input_handle, target.open("xb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=8 * 1024 * 1024)
        output_handle.flush()
        os.fsync(output_handle.fileno())


def _validate_model_config(
    path: Path, manifest: Mapping[str, Any], model_id: str
) -> None:
    contracts = manifest.get("architecture_contract")
    if not isinstance(contracts, Mapping) or not isinstance(
        contracts.get(model_id), Mapping
    ):
        raise CheckpointPreflightError(
            f"checkpoint manifest lacks the {model_id} architecture contract"
        )
    expected = contracts[model_id]
    required_fields = {
        "architectures",
        "hidden_size",
        "intermediate_size",
        "max_position_embeddings",
        "model_type",
        "num_attention_heads",
        "num_hidden_layers",
        "pad_token_id",
        "vocab_size",
    }
    if set(expected) != required_fields:
        raise CheckpointPreflightError(
            f"{model_id} architecture contract fields differ"
        )
    try:
        observed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CheckpointPreflightError(
            f"{model_id} config is invalid JSON: {error}"
        ) from error
    if not isinstance(observed, Mapping) or any(
        observed.get(field) != expected[field] for field in required_fields
    ):
        raise CheckpointPreflightError(
            f"{model_id} config differs from the frozen architecture"
        )
    positive_integer_fields = required_fields.difference(
        {"architectures", "model_type", "pad_token_id"}
    )
    if any(
        isinstance(expected[field], bool)
        or not isinstance(expected[field], int)
        or expected[field] <= 0
        for field in positive_integer_fields
    ):
        raise CheckpointPreflightError(
            f"{model_id} architecture dimensions are invalid"
        )


def _sanitize_dictionary(
    source: Path, target: Path, spec: Mapping[str, Any]
) -> dict[str, Any]:
    before = source.stat(follow_symlinks=False)
    payload = source.read_bytes()
    after = source.stat(follow_symlinks=False)
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ) or source.is_symlink():
        raise CheckpointPreflightError(
            f"checkpoint dictionary changed while reading: {spec['path']}"
        )
    if len(payload) != spec["size_bytes"] or sha256_bytes(payload) != spec["sha256"]:
        raise CheckpointPreflightError(
            f"checkpoint dictionary content differs: {spec['path']}"
        )
    mapping = load_restricted_plain_mapping(payload)
    target.parent.mkdir(parents=True, exist_ok=True)
    write_bytes_exclusive(
        target,
        pickle.dumps(mapping, protocol=4, fix_imports=False),
    )
    json_target = target.with_suffix(target.suffix + ".json")
    write_json_exclusive(json_target, mapping)
    return {
        "entry_count": len(mapping),
        "canonical_mapping_sha256": canonical_sha256(mapping),
        "sanitized_pickle_path": target.as_posix(),
        "canonical_json_path": json_target.as_posix(),
    }


def stage_geneformer_bundle(
    *,
    model_id: str,
    source_root: str | Path,
    output_root: str | Path,
    manifest_path: str | Path,
) -> Path:
    """Verify and publish one read-only, execution-safe Geneformer bundle."""

    manifest_file = reject_symlink_components(
        Path(manifest_path), label="checkpoint manifest"
    )
    manifest = _load_manifest(manifest_file)
    required = _required_artifacts(manifest, model_id)
    source = reject_symlink_components(Path(source_root), label="checkpoint source root")
    destination_root = reject_symlink_components(
        Path(output_root), label="checkpoint output root"
    )
    if not source.is_dir() or not destination_root.is_dir():
        raise CheckpointPreflightError(
            "checkpoint source and output roots must already exist as directories"
        )
    manifest_sha256 = sha256_file(manifest_file)
    target = destination_root / f"{model_id}--{manifest_sha256[:16]}"
    if target.exists() or target.is_symlink():
        raise CheckpointPreflightError(
            f"refusing to replace existing checkpoint bundle: {target}"
        )

    temporary = Path(tempfile.mkdtemp(prefix=f".{model_id}.", dir=destination_root))
    try:
        staged: list[dict[str, Any]] = []
        dictionaries: list[dict[str, Any]] = []
        for spec in required:
            source_path = _verified_source(source / spec["path"], spec)
            staged_path = temporary / "upstream" / spec["path"]
            if spec["path"].endswith(".pkl"):
                dictionary = _sanitize_dictionary(
                    source_path,
                    temporary / "sanitized" / spec["path"],
                    spec,
                )
                dictionaries.append(
                    {
                        "source_path": spec["path"],
                        **{
                            key: (
                                Path(value).relative_to(temporary).as_posix()
                                if key.endswith("_path")
                                else value
                            )
                            for key, value in dictionary.items()
                        },
                    }
                )
            else:
                _copy_exclusive(source_path, staged_path)
                if (
                    staged_path.stat().st_size != spec["size_bytes"]
                    or sha256_file(staged_path) != spec["sha256"]
                ):
                    raise CheckpointPreflightError(
                        f"staged artifact differs after copying: {spec['path']}"
                    )
                if Path(spec["path"]).name == "config.json":
                    _validate_model_config(staged_path, manifest, model_id)
            staged.append(dict(spec))

        bundle = {
            "schema_version": "masld-bench-geneformer-bundle-v1",
            "model_id": model_id,
            "repository": manifest["repository"],
            "repository_revision": manifest["repository_revision"],
            "source_manifest_sha256": manifest_sha256,
            "required_upstream_artifacts": staged,
            "sanitized_dictionaries": dictionaries,
            "weight_content_executed": False,
            "pickle_global_or_class_construction_allowed": False,
            "downloads_performed": False,
        }
        write_json_exclusive(temporary / "bundle.json", bundle)
        freeze_tree(
            temporary,
            metadata={
                "artifact_class": "geneformer_checkpoint_bundle",
                "model_id": model_id,
                "source_manifest_sha256": manifest_sha256,
            },
        )
        verify_frozen_tree(temporary)
        publish_directory_noreplace(temporary, target)
        verify_frozen_tree(target)
    except (ArtifactError, CheckpointPreflightError, OSError, ValueError) as error:
        shutil.rmtree(temporary, ignore_errors=True)
        if isinstance(error, CheckpointPreflightError):
            raise
        raise CheckpointPreflightError(f"Geneformer staging failed: {error}") from error
    return target


__all__ = [
    "CheckpointPreflightError",
    "load_restricted_plain_mapping",
    "stage_geneformer_bundle",
]
