"""Strict, immutable contracts for the MASLD model benchmark.

The registry is an executable scientific contract rather than a permissive
configuration bag.  Every loader rejects unknown keys, every enum is explicit,
and file-bearing records verify size and SHA-256 before use.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date, datetime
from enum import Enum
import json
from pathlib import Path
import re
import tomllib
from types import MappingProxyType
from typing import Any, ClassVar, Mapping, TypeVar

from .artifacts import ArtifactError, reject_symlink_components
from .hashing import (
    HashingError,
    canonical_json,
    canonical_sha256,
    canonicalize,
    require_sha256,
    sha256_file,
)


class ContractError(ValueError):
    """Raised when a benchmark contract is incomplete or ambiguous."""


class PairingState(str, Enum):
    SAME_MOLECULE = "same_molecule"
    SAME_CELL = "same_cell"
    SAME_NUCLEUS = "same_nucleus"
    SAME_SECTION = "same_section"
    ADJACENT_SECTION = "adjacent_section"
    SAME_SAMPLE_DIFFERENT_ALIQUOT = "same_sample_different_aliquot"
    SAME_DONOR_DIFFERENT_TISSUE = "same_donor_different_tissue"
    SAME_STUDY_UNPAIRED = "same_study_unpaired"


class MissingState(str, Enum):
    OBSERVED = "observed"
    STRUCTURALLY_MISSING = "structurally_missing"
    NOT_APPLICABLE = "not_applicable"
    BELOW_QC = "below_qc"
    UNAVAILABLE_PERMISSION = "unavailable_permission"
    JOIN_UNRESOLVED = "join_unresolved"
    WITHHELD_SEALED = "withheld_sealed"
    DERIVABLE_NOT_PROCESSED = "derivable_not_processed"


class ExposureState(str, Enum):
    CLEAN_DECLARED = "clean_declared"
    TARGET_LABEL_UNEXPOSED = "target_label_unexposed"
    ENCODER_SEEN = "encoder_seen"
    CONTINUAL_SEEN = "continual_seen"
    REFERENCE_ONLY = "reference_only"
    DOWNSTREAM_DEMO = "downstream_demo"
    UNKNOWN = "unknown"


class AccessState(str, Enum):
    PUBLIC = "public"
    CONTROLLED = "controlled"
    CONDITIONAL = "conditional"


class ReleaseClass(str, Enum):
    CANDIDATE = "candidate"
    RESTRICTED_COMPARATOR = "restricted_comparator"
    BLOCKED_TERMS = "blocked_terms"
    DEFERRED = "deferred"


class ReleaseState(str, Enum):
    WORKING = "working"
    CANDIDATE = "candidate"
    RESTRICTED_COMPARATOR = "restricted_comparator"
    BLOCKED_TERMS = "blocked_terms"
    DEFERRED = "deferred"
    SELECTED = "selected"
    PROMOTED = "promoted"
    RELEASED = "released"
    SUPERSEDED = "superseded"
    REJECTED = "rejected"
    BLOCKED = "blocked"


class RunState(str, Enum):
    PLANNED = "planned"
    SUBMITTED = "submitted"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    INVALID_OUTPUT = "invalid_output"
    FAILED_SOFTWARE = "failed_software"
    FAILED_RESOURCE = "failed_resource"
    BLOCKED_UPSTREAM = "blocked_upstream"
    BLOCKED_LICENSE = "blocked_license"
    BLOCKED_DATA = "blocked_data"
    SKIPPED_GATE = "skipped_gate"


class DatasetRole(str, Enum):
    TRAIN_DEVELOPMENT = "train_development"
    EXTERNAL_DEVELOPMENT = "external_development"
    WITHHELD_SEALED = "withheld_sealed"
    MECHANISTIC_DIAGNOSTIC = "mechanistic_diagnostic"
    TECHNICAL_PILOT = "technical_pilot"
    CROSS_SPECIES_STRESS = "cross_species_stress"
    BLOCKED = "blocked"
    DEFERRED = "deferred"


class DatasetStatus(str, Enum):
    AVAILABLE = "available"
    WITHHELD_SEALED = "withheld_sealed"
    BLOCKED = "blocked"
    DEFERRED = "deferred"


class TaskStatus(str, Enum):
    CANDIDATE = "candidate"
    BLOCKED = "blocked"
    DEFERRED = "deferred"


class ImplementationType(str, Enum):
    PRETRAINED_ADAPTER = "pretrained_adapter"
    TRAIN_LOCALLY = "train_locally"
    CLASSICAL_BASELINE = "classical_baseline"
    GRAPH_BASELINE = "graph_baseline"


_ContractT = TypeVar("_ContractT", bound="StrictContract")
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_IDENTIFIER_PATTERN = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")


def _contains_unresolved_marker(value: str) -> bool:
    """Recognize exact sentinels and delimited embedded markers."""

    stripped = value.strip()
    if stripped.upper() in {"", "UNKNOWN", "UNRESOLVED", "TBD", "LATEST"}:
        return True
    return bool(
        re.search(
            r"(?:^|[^A-Z0-9])(?:UNKNOWN|UNRESOLVED|TBD|LATEST)(?:$|[^A-Z0-9])",
            stripped.upper(),
        )
    )


def _strict_mapping(
    value: object,
    *,
    allowed: set[str] | frozenset[str],
    required: set[str] | frozenset[str],
    label: str,
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ContractError(f"{label} must be an object/table")
    if any(not isinstance(key, str) for key in value):
        raise ContractError(f"{label} keys must be strings")
    unknown = sorted(set(value) - set(allowed))
    if unknown:
        raise ContractError(f"{label} has unknown field(s): {', '.join(unknown)}")
    missing = sorted(set(required) - set(value))
    if missing:
        raise ContractError(f"{label} is missing field(s): {', '.join(missing)}")
    return dict(value)


def _as_string(value: object, field_name: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise ContractError(f"{field_name} must be a string")
    if not allow_empty and not value.strip():
        raise ContractError(f"{field_name} must not be empty")
    if "\x00" in value:
        raise ContractError(f"{field_name} must not contain NUL")
    return value


def _optional_string(value: object, field_name: str) -> str | None:
    if value is None:
        return None
    return _as_string(value, field_name)


def _as_bool(value: object, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise ContractError(f"{field_name} must be boolean")
    return value


def _as_int(value: object, field_name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractError(f"{field_name} must be an integer")
    if value < minimum:
        raise ContractError(f"{field_name} must be >= {minimum}")
    return value


def _optional_int(
    value: object, field_name: str, *, minimum: int = 0
) -> int | None:
    if value is None:
        return None
    return _as_int(value, field_name, minimum=minimum)


def _optional_int_or_unresolved(
    value: object, field_name: str, *, minimum: int = 0
) -> int | str | None:
    if value == "UNRESOLVED":
        return "UNRESOLVED"
    return _optional_int(value, field_name, minimum=minimum)


def _as_enum(value: object, enum_type: type[Enum], field_name: str) -> Any:
    if isinstance(value, enum_type):
        return value
    if not isinstance(value, str):
        raise ContractError(f"{field_name} must be a string")
    try:
        return enum_type(value)
    except ValueError as error:
        choices = ", ".join(repr(item.value) for item in enum_type)
        raise ContractError(f"{field_name} must be one of: {choices}") from error


def _as_string_tuple(
    value: object,
    field_name: str,
    *,
    allow_empty: bool = True,
) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        raise ContractError(f"{field_name} must be an array")
    result = tuple(
        _as_string(item, f"{field_name}[{index}]")
        for index, item in enumerate(value)
    )
    if not allow_empty and not result:
        raise ContractError(f"{field_name} must not be empty")
    if len(set(result)) != len(result):
        raise ContractError(f"{field_name} must not contain duplicates")
    return result


def _as_enum_tuple(
    value: object,
    enum_type: type[Enum],
    field_name: str,
    *,
    allow_empty: bool = True,
) -> tuple[Any, ...]:
    if not isinstance(value, (list, tuple)):
        raise ContractError(f"{field_name} must be an array")
    result = tuple(
        _as_enum(item, enum_type, f"{field_name}[{index}]")
        for index, item in enumerate(value)
    )
    if not allow_empty and not result:
        raise ContractError(f"{field_name} must not be empty")
    if len(set(result)) != len(result):
        raise ContractError(f"{field_name} must not contain duplicates")
    return result


def _freeze_json(value: object, field_name: str) -> Any:
    try:
        normalized = canonicalize(value)
    except (HashingError, TypeError, ValueError) as error:
        raise ContractError(f"{field_name} must contain only JSON values: {error}") from error

    def freeze(item: Any) -> Any:
        if isinstance(item, dict):
            return MappingProxyType({key: freeze(child) for key, child in item.items()})
        if isinstance(item, list):
            return tuple(freeze(child) for child in item)
        return item

    return freeze(normalized)


def _as_metadata(value: object, field_name: str = "metadata") -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ContractError(f"{field_name} must be an object/table")
    return _freeze_json(value, field_name)


def _as_sha256(value: object, field_name: str) -> str:
    try:
        return require_sha256(value, field_name=field_name)
    except HashingError as error:
        raise ContractError(str(error)) from error


def _as_sha256_or_unresolved(value: object, field_name: str) -> str:
    if value in {"UNRESOLVED", "NOT_APPLICABLE"}:
        return str(value)
    return _as_sha256(value, field_name)


def _as_verified_date(value: object, field_name: str) -> str:
    result = _as_string(value, field_name)
    try:
        date.fromisoformat(result)
    except ValueError as error:
        raise ContractError(f"{field_name} must be an ISO 8601 date (YYYY-MM-DD)") from error
    return result


def _as_utc_timestamp(value: object, field_name: str) -> str:
    result = _as_string(value, field_name)
    if not result.endswith("Z"):
        raise ContractError(f"{field_name} must end in Z")
    try:
        datetime.fromisoformat(result[:-1] + "+00:00")
    except ValueError as error:
        raise ContractError(f"{field_name} must be an ISO 8601 UTC timestamp") from error
    return result


def _set_frozen(instance: object, field_name: str, value: object) -> None:
    object.__setattr__(instance, field_name, value)


class StrictContract:
    """Serialization helpers shared by all strict dataclass contracts."""

    @classmethod
    def from_dict(cls: type[_ContractT], value: Mapping[str, Any]) -> _ContractT:
        raise NotImplementedError

    @classmethod
    def from_json(cls: type[_ContractT], document: str | bytes) -> _ContractT:
        try:
            value = json.loads(document)
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            raise ContractError(f"invalid JSON for {cls.__name__}: {error}") from error
        return cls.from_dict(value)

    @classmethod
    def from_toml(cls: type[_ContractT], document: str) -> _ContractT:
        try:
            value = tomllib.loads(document)
        except tomllib.TOMLDecodeError as error:
            raise ContractError(f"invalid TOML for {cls.__name__}: {error}") from error
        return cls.from_dict(value)

    @classmethod
    def load_json(cls: type[_ContractT], path: str | Path) -> _ContractT:
        source = Path(path)
        try:
            document = source.read_bytes()
        except OSError as error:
            raise ContractError(f"cannot read {source}: {error}") from error
        return cls.from_json(document)

    @classmethod
    def load_toml(cls: type[_ContractT], path: str | Path) -> _ContractT:
        source = Path(path)
        try:
            document = source.read_text(encoding="utf-8")
        except OSError as error:
            raise ContractError(f"cannot read {source}: {error}") from error
        return cls.from_toml(document)

    def to_dict(self) -> dict[str, Any]:
        result = canonicalize(self)
        if not isinstance(result, dict):
            raise ContractError("contract did not serialize to a JSON object")
        return result

    def as_dict(self) -> dict[str, Any]:
        return self.to_dict()

    def to_json(self) -> str:
        return canonical_json(self.to_dict())


@dataclass(frozen=True, slots=True)
class ArtifactRef(StrictContract):
    path: str
    sha256: str
    size_bytes: int
    media_type: str | None = None
    role: str | None = None

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {"path", "sha256", "size_bytes", "media_type", "role"}
    )
    _REQUIRED: ClassVar[frozenset[str]] = frozenset(
        {"path", "sha256", "size_bytes"}
    )

    def __post_init__(self) -> None:
        _set_frozen(self, "path", _as_string(self.path, "ArtifactRef.path"))
        _set_frozen(self, "sha256", _as_sha256(self.sha256, "ArtifactRef.sha256"))
        _set_frozen(
            self,
            "size_bytes",
            _as_int(self.size_bytes, "ArtifactRef.size_bytes", minimum=0),
        )
        _set_frozen(
            self,
            "media_type",
            _optional_string(self.media_type, "ArtifactRef.media_type"),
        )
        _set_frozen(self, "role", _optional_string(self.role, "ArtifactRef.role"))

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ArtifactRef":
        raw = _strict_mapping(
            value, allowed=cls._FIELDS, required=cls._REQUIRED, label="ArtifactRef"
        )
        return cls(
            path=raw["path"],
            sha256=raw["sha256"],
            size_bytes=raw["size_bytes"],
            media_type=raw.get("media_type"),
            role=raw.get("role"),
        )

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        *,
        relative_to: str | Path | None = None,
        media_type: str | None = None,
        role: str | None = None,
    ) -> "ArtifactRef":
        try:
            source = reject_symlink_components(Path(path), label="artifact path")
        except ArtifactError as error:
            raise ContractError(str(error)) from error
        if not source.is_file():
            raise ContractError(f"artifact is not a regular file: {source}")
        display = source
        if relative_to is not None:
            try:
                relative_root = reject_symlink_components(
                    Path(relative_to), label="artifact root"
                )
            except ArtifactError as error:
                raise ContractError(str(error)) from error
            try:
                display = source.resolve().relative_to(relative_root.resolve())
            except ValueError as error:
                raise ContractError(
                    f"artifact {source} is outside root {relative_to}"
                ) from error
        return cls(
            path=display.as_posix(),
            sha256=sha256_file(source),
            size_bytes=source.stat().st_size,
            media_type=media_type,
            role=role,
        )

    def resolve_path(
        self,
        root: str | Path | None = None,
        *,
        require_relative: bool = False,
    ) -> Path:
        configured = Path(self.path)
        if require_relative and configured.is_absolute():
            raise ContractError(f"artifact path must be relative: {self.path}")
        try:
            configured_root = (
                reject_symlink_components(Path(root), label="artifact root")
                if root is not None
                else None
            )
        except ArtifactError as error:
            raise ContractError(str(error)) from error
        base = configured_root.resolve() if configured_root is not None else None
        candidate = (
            configured
            if configured.is_absolute()
            else (base / configured if base else configured)
        )
        try:
            candidate = reject_symlink_components(candidate, label="artifact path")
        except ArtifactError as error:
            raise ContractError(str(error)) from error
        try:
            resolved = candidate.resolve(strict=True)
        except OSError as error:
            raise ContractError(f"artifact is missing or unreadable: {candidate}") from error
        if base is not None:
            try:
                resolved.relative_to(base)
            except ValueError as error:
                raise ContractError(
                    f"artifact escapes declared root {base}: {self.path}"
                ) from error
        if not resolved.is_file():
            raise ContractError(f"artifact is not a regular file: {resolved}")
        return resolved

    def validate(
        self,
        root: str | Path | None = None,
        *,
        require_relative: bool = False,
    ) -> Path:
        source = self.resolve_path(root, require_relative=require_relative)
        before = source.stat(follow_symlinks=False)
        actual_sha256 = sha256_file(source)
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
        ):
            raise ContractError(f"artifact changed while hashing: {source}")
        actual_size = after.st_size
        if actual_size != self.size_bytes:
            raise ContractError(
                f"artifact size mismatch for {source}: expected {self.size_bytes}, got {actual_size}"
            )
        if actual_sha256 != self.sha256:
            raise ContractError(
                f"artifact SHA-256 mismatch for {source}: expected {self.sha256}, got {actual_sha256}"
            )
        return source


@dataclass(frozen=True, slots=True)
class DatasetProvenance(StrictContract):
    source_url: str
    verified_date: str
    local_exposure: str

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {"source_url", "verified_date", "local_exposure"}
    )

    def __post_init__(self) -> None:
        _set_frozen(
            self, "source_url", _as_string(self.source_url, "provenance.source_url")
        )
        _set_frozen(
            self,
            "verified_date",
            _as_verified_date(self.verified_date, "provenance.verified_date"),
        )
        _set_frozen(
            self,
            "local_exposure",
            _as_string(self.local_exposure, "provenance.local_exposure"),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DatasetProvenance":
        raw = _strict_mapping(
            value, allowed=cls._FIELDS, required=cls._FIELDS, label="DatasetProvenance"
        )
        return cls(**raw)


@dataclass(frozen=True, slots=True)
class DatasetSplit(StrictContract):
    group_key: str
    outer_role: str
    label_visibility: str

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {"group_key", "outer_role", "label_visibility"}
    )

    def __post_init__(self) -> None:
        for field_name in self._FIELDS:
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"split.{field_name}"),
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DatasetSplit":
        raw = _strict_mapping(
            value, allowed=cls._FIELDS, required=cls._FIELDS, label="DatasetSplit"
        )
        return cls(**raw)


@dataclass(frozen=True, slots=True)
class DatasetActivationContract(StrictContract):
    """Checksummed admission evidence for one immutable dataset version.

    Registry records may describe a dataset before it is executable.  A dataset
    becomes runnable only after this contract binds the rights, topology, donor
    join, labels, QC, reference, and concrete artifact inventory reviewed for
    that exact version.  Project-sealed datasets keep this contract unresolved
    until the isolated evaluator activates anonymized features.
    """

    schema_version: str = "masld-bench-dataset-activation-v1"
    dataset_version: str = "UNRESOLVED"
    rights_sha256: str = "UNRESOLVED"
    topology_sha256: str = "UNRESOLVED"
    donor_join_sha256: str = "UNRESOLVED"
    labels_sha256: str = "UNRESOLVED"
    qc_sha256: str = "UNRESOLVED"
    reference_sha256: str = "UNRESOLVED"
    exposure_audit_sha256: str = "UNRESOLVED"
    evidence_artifacts: tuple[ArtifactRef, ...] = ()
    artifact_manifest: ArtifactRef | None = None
    ready: bool = False

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "dataset_version",
            "rights_sha256",
            "topology_sha256",
            "donor_join_sha256",
            "labels_sha256",
            "qc_sha256",
            "reference_sha256",
            "exposure_audit_sha256",
            "evidence_artifacts",
            "artifact_manifest",
            "ready",
        }
    )
    _AUTHORITY_FIELDS: ClassVar[tuple[str, ...]] = (
        "rights",
        "topology",
        "donor_join",
        "labels",
        "qc",
        "reference",
        "exposure_audit",
    )
    def __post_init__(self) -> None:
        _set_frozen(
            self,
            "schema_version",
            _as_string(self.schema_version, "DatasetActivationContract.schema_version"),
        )
        if self.schema_version != "masld-bench-dataset-activation-v1":
            raise ContractError(
                "DatasetActivationContract.schema_version must be "
                "masld-bench-dataset-activation-v1"
            )
        _set_frozen(
            self,
            "dataset_version",
            _as_string(self.dataset_version, "DatasetActivationContract.dataset_version"),
        )
        for field_name in (
            "rights_sha256",
            "topology_sha256",
            "donor_join_sha256",
            "labels_sha256",
            "qc_sha256",
            "reference_sha256",
            "exposure_audit_sha256",
        ):
            _set_frozen(
                self,
                field_name,
                _as_sha256_or_unresolved(
                    getattr(self, field_name),
                    f"DatasetActivationContract.{field_name}",
                ),
            )
        if self.artifact_manifest is not None and not isinstance(
            self.artifact_manifest, ArtifactRef
        ):
            raise ContractError(
                "DatasetActivationContract.artifact_manifest must be ArtifactRef or null"
            )
        if not isinstance(self.evidence_artifacts, (tuple, list)):
            raise ContractError(
                "DatasetActivationContract.evidence_artifacts must be an array"
            )
        evidence_artifacts = tuple(self.evidence_artifacts)
        if any(not isinstance(item, ArtifactRef) for item in evidence_artifacts):
            raise ContractError(
                "DatasetActivationContract.evidence_artifacts entries must be ArtifactRef"
            )
        evidence_paths = [item.path for item in evidence_artifacts]
        evidence_roles = [item.role for item in evidence_artifacts]
        if len(set(evidence_paths)) != len(evidence_paths):
            raise ContractError(
                "DatasetActivationContract.evidence_artifacts paths must be unique"
            )
        if any(role is None for role in evidence_roles) or len(set(evidence_roles)) != len(
            evidence_roles
        ):
            raise ContractError(
                "DatasetActivationContract.evidence_artifacts require unique roles"
            )
        suffix_to_artifact: dict[str, ArtifactRef] = {}
        for artifact in evidence_artifacts:
            assert artifact.role is not None
            parts = artifact.role.split(":")
            if (
                len(parts) != 3
                or parts[0] != "dataset_authority"
                or not parts[1]
                or parts[2] not in self._AUTHORITY_FIELDS
            ):
                raise ContractError(
                    "DatasetActivationContract evidence role must be "
                    "dataset_authority:<dataset_id>:<authority>"
                )
            suffix_to_artifact[parts[2]] = artifact
        if evidence_artifacts and set(suffix_to_artifact) != set(self._AUTHORITY_FIELDS):
            raise ContractError(
                "DatasetActivationContract.evidence_artifacts must cover every and only "
                "rights, topology, donor_join, labels, qc, reference, and exposure_audit"
            )
        for authority, artifact in suffix_to_artifact.items():
            if artifact.sha256 != getattr(self, f"{authority}_sha256"):
                raise ContractError(
                    "DatasetActivationContract authority artifact SHA-256 differs from "
                    f"{authority}_sha256"
                )
        _set_frozen(self, "evidence_artifacts", evidence_artifacts)
        _set_frozen(
            self,
            "ready",
            _as_bool(self.ready, "DatasetActivationContract.ready"),
        )
        unresolved = (
            self.dataset_version == "UNRESOLVED"
            or any(
                getattr(self, field_name) == "UNRESOLVED"
                for field_name in (
                    "rights_sha256",
                    "topology_sha256",
                    "donor_join_sha256",
                    "labels_sha256",
                    "qc_sha256",
                    "reference_sha256",
                    "exposure_audit_sha256",
                )
            )
            or self.artifact_manifest is None
            or len(evidence_artifacts) != len(self._AUTHORITY_FIELDS)
        )
        if self.ready and unresolved:
            raise ContractError(
                "ready dataset activation requires a version, every authority hash, "
                "and an artifact manifest"
            )
        if not self.ready and not unresolved:
            raise ContractError(
                "complete dataset activation evidence must set ready=true"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DatasetActivationContract":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS,
            label="DatasetActivationContract",
        )
        artifact = raw["artifact_manifest"]
        evidence = raw["evidence_artifacts"]
        return cls(
            schema_version=raw["schema_version"],
            dataset_version=raw["dataset_version"],
            rights_sha256=raw["rights_sha256"],
            topology_sha256=raw["topology_sha256"],
            donor_join_sha256=raw["donor_join_sha256"],
            labels_sha256=raw["labels_sha256"],
            qc_sha256=raw["qc_sha256"],
            reference_sha256=raw["reference_sha256"],
            exposure_audit_sha256=raw["exposure_audit_sha256"],
            evidence_artifacts=tuple(
                ArtifactRef.from_dict(item) for item in evidence
            ),
            artifact_manifest=(
                ArtifactRef.from_dict(artifact) if artifact is not None else None
            ),
            ready=raw["ready"],
        )


@dataclass(frozen=True, slots=True)
class DatasetViewContract(StrictContract):
    """Immutable, purpose-limited view of one registered parent dataset.

    A view can make a small compatibility fixture executable without claiming
    that the complete parent dataset has passed activation.  The permission is
    intentionally narrow: v1 views are development-visible smoke inputs only,
    and carry the same seven scientific authorities as a full activation.
    """

    schema_version: str
    view_id: str
    parent_dataset_id: str
    parent_registry_sha256: str
    purpose: str
    allowed_waves: tuple[str, ...]
    row_count: int
    biological_unit_count: int
    modalities: tuple[str, ...]
    pairing_levels: tuple[PairingState, ...]
    label_visibility: str
    selection_policy: str
    selection_seed: int
    selection_outcomes_used: bool
    sealed_outcomes_used: bool
    class_counts: Mapping[str, Any]
    artifact_manifest: ArtifactRef
    complete_receipt: ArtifactRef
    data_artifact: ArtifactRef
    selection_artifact: ArtifactRef
    ontology_artifact: ArtifactRef
    derivation_artifact: ArtifactRef
    authority_artifacts: tuple[ArtifactRef, ...]
    ready: bool

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "view_id",
            "parent_dataset_id",
            "parent_registry_sha256",
            "purpose",
            "allowed_waves",
            "row_count",
            "biological_unit_count",
            "modalities",
            "pairing_levels",
            "label_visibility",
            "selection_policy",
            "selection_seed",
            "selection_outcomes_used",
            "sealed_outcomes_used",
            "class_counts",
            "artifact_manifest",
            "complete_receipt",
            "data_artifact",
            "selection_artifact",
            "ontology_artifact",
            "derivation_artifact",
            "authority_artifacts",
            "ready",
        }
    )
    _AUTHORITY_FIELDS: ClassVar[tuple[str, ...]] = (
        "rights",
        "topology",
        "donor_join",
        "labels",
        "qc",
        "reference",
        "exposure_audit",
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "view_id",
            "parent_dataset_id",
            "purpose",
            "label_visibility",
            "selection_policy",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(
                    getattr(self, field_name), f"DatasetViewContract.{field_name}"
                ),
            )
        if self.schema_version != "masld-bench-dataset-view-v1":
            raise ContractError(
                "DatasetViewContract.schema_version must be "
                "masld-bench-dataset-view-v1"
            )
        if self.purpose != "compatibility_smoke":
            raise ContractError(
                "DatasetViewContract.purpose must be compatibility_smoke"
            )
        _set_frozen(
            self,
            "parent_registry_sha256",
            _as_sha256(
                self.parent_registry_sha256,
                "DatasetViewContract.parent_registry_sha256",
            ),
        )
        _set_frozen(
            self,
            "allowed_waves",
            _as_string_tuple(
                self.allowed_waves,
                "DatasetViewContract.allowed_waves",
                allow_empty=False,
            ),
        )
        if self.allowed_waves != ("smoke",):
            raise ContractError(
                "DatasetViewContract.allowed_waves must be exactly ['smoke']"
            )
        for field_name in ("row_count", "biological_unit_count"):
            _set_frozen(
                self,
                field_name,
                _as_int(
                    getattr(self, field_name),
                    f"DatasetViewContract.{field_name}",
                    minimum=1,
                ),
            )
        if self.biological_unit_count > self.row_count:
            raise ContractError(
                "DatasetViewContract biological_unit_count cannot exceed row_count"
            )
        _set_frozen(
            self,
            "modalities",
            _as_string_tuple(
                self.modalities, "DatasetViewContract.modalities", allow_empty=False
            ),
        )
        _set_frozen(
            self,
            "pairing_levels",
            _as_enum_tuple(
                self.pairing_levels,
                PairingState,
                "DatasetViewContract.pairing_levels",
                allow_empty=False,
            ),
        )
        if self.label_visibility != "development_visible":
            raise ContractError(
                "DatasetViewContract must keep labels development_visible"
            )
        _set_frozen(
            self,
            "selection_seed",
            _as_int(
                self.selection_seed,
                "DatasetViewContract.selection_seed",
                minimum=0,
            ),
        )
        for field_name in ("selection_outcomes_used", "sealed_outcomes_used", "ready"):
            _set_frozen(
                self,
                field_name,
                _as_bool(
                    getattr(self, field_name), f"DatasetViewContract.{field_name}"
                ),
            )
        if self.selection_outcomes_used or self.sealed_outcomes_used:
            raise ContractError(
                "DatasetViewContract may not use selection or sealed outcomes"
            )
        _set_frozen(
            self,
            "class_counts",
            _as_metadata(self.class_counts, "DatasetViewContract.class_counts"),
        )
        counts = dict(self.class_counts)
        if (
            not counts
            or any(
                not isinstance(label, str)
                or not label
                or isinstance(count, bool)
                or not isinstance(count, int)
                or count < 1
                for label, count in counts.items()
            )
            or sum(counts.values()) != self.row_count
        ):
            raise ContractError(
                "DatasetViewContract.class_counts must be positive and sum to row_count"
            )
        artifact_roles = {
            "artifact_manifest": f"dataset_view_manifest:{self.view_id}",
            "complete_receipt": f"dataset_view_complete:{self.view_id}",
            "data_artifact": f"dataset_view_data:{self.view_id}",
            "selection_artifact": f"dataset_view_selection:{self.view_id}",
            "ontology_artifact": f"dataset_view_ontology:{self.view_id}",
            "derivation_artifact": f"dataset_view_derivation:{self.view_id}",
        }
        direct_artifacts: list[ArtifactRef] = []
        for field_name, expected_role in artifact_roles.items():
            artifact = getattr(self, field_name)
            if not isinstance(artifact, ArtifactRef):
                raise ContractError(
                    f"DatasetViewContract.{field_name} must be ArtifactRef"
                )
            if artifact.role != expected_role:
                raise ContractError(
                    f"DatasetViewContract.{field_name} must use role {expected_role}"
                )
            direct_artifacts.append(artifact)
        if not isinstance(self.authority_artifacts, (tuple, list)):
            raise ContractError(
                "DatasetViewContract.authority_artifacts must be an array"
            )
        authorities = tuple(self.authority_artifacts)
        if any(not isinstance(item, ArtifactRef) for item in authorities):
            raise ContractError(
                "DatasetViewContract.authority_artifacts entries must be ArtifactRef"
            )
        expected_authority_roles = {
            f"dataset_view_authority:{self.view_id}:{authority}"
            for authority in self._AUTHORITY_FIELDS
        }
        if {item.role for item in authorities} != expected_authority_roles:
            raise ContractError(
                "DatasetViewContract authority_artifacts must cover every and only "
                "rights, topology, donor_join, labels, qc, reference, and exposure_audit"
            )
        all_paths = [item.path for item in (*direct_artifacts, *authorities)]
        if len(set(all_paths)) != len(all_paths):
            raise ContractError("DatasetViewContract artifact paths must be unique")
        _set_frozen(self, "authority_artifacts", authorities)
        if not self.ready:
            raise ContractError(
                "registered DatasetViewContract must be complete and ready"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DatasetViewContract":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS,
            label="DatasetViewContract",
        )
        return cls(
            schema_version=raw["schema_version"],
            view_id=raw["view_id"],
            parent_dataset_id=raw["parent_dataset_id"],
            parent_registry_sha256=raw["parent_registry_sha256"],
            purpose=raw["purpose"],
            allowed_waves=raw["allowed_waves"],
            row_count=raw["row_count"],
            biological_unit_count=raw["biological_unit_count"],
            modalities=raw["modalities"],
            pairing_levels=raw["pairing_levels"],
            label_visibility=raw["label_visibility"],
            selection_policy=raw["selection_policy"],
            selection_seed=raw["selection_seed"],
            selection_outcomes_used=raw["selection_outcomes_used"],
            sealed_outcomes_used=raw["sealed_outcomes_used"],
            class_counts=raw["class_counts"],
            artifact_manifest=ArtifactRef.from_dict(raw["artifact_manifest"]),
            complete_receipt=ArtifactRef.from_dict(raw["complete_receipt"]),
            data_artifact=ArtifactRef.from_dict(raw["data_artifact"]),
            selection_artifact=ArtifactRef.from_dict(raw["selection_artifact"]),
            ontology_artifact=ArtifactRef.from_dict(raw["ontology_artifact"]),
            derivation_artifact=ArtifactRef.from_dict(raw["derivation_artifact"]),
            authority_artifacts=tuple(
                ArtifactRef.from_dict(item) for item in raw["authority_artifacts"]
            ),
            ready=raw["ready"],
        )

    @property
    def supporting_artifacts(self) -> tuple[ArtifactRef, ...]:
        return (
            self.complete_receipt,
            self.data_artifact,
            self.selection_artifact,
            self.ontology_artifact,
            self.derivation_artifact,
            *self.authority_artifacts,
        )


@dataclass(frozen=True, slots=True)
class PredictionFirstOutcomeContract(StrictContract):
    """Policy for visible design features with outcomes withheld until commit."""

    schema_version: str
    policy_id: str
    scope: str
    design_visibility: str
    outcome_visibility: str
    prediction_scope: str
    unit_of_replication: str
    prediction_universe_sha256: str
    prediction_universe_role: str
    excluded_endpoint_ids: tuple[str, ...]
    forbid_fit: bool
    forbid_model_selection: bool
    require_prediction_commit: bool

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "policy_id",
            "scope",
            "design_visibility",
            "outcome_visibility",
            "prediction_scope",
            "unit_of_replication",
            "prediction_universe_sha256",
            "prediction_universe_role",
            "excluded_endpoint_ids",
            "forbid_fit",
            "forbid_model_selection",
            "require_prediction_commit",
        }
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "policy_id",
            "scope",
            "design_visibility",
            "outcome_visibility",
            "prediction_scope",
            "unit_of_replication",
            "prediction_universe_role",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(
                    getattr(self, field_name),
                    f"PredictionFirstOutcomeContract.{field_name}",
                ),
            )
        _set_frozen(
            self,
            "prediction_universe_sha256",
            _as_sha256_or_unresolved(
                self.prediction_universe_sha256,
                "PredictionFirstOutcomeContract.prediction_universe_sha256",
            ),
        )
        _set_frozen(
            self,
            "excluded_endpoint_ids",
            _as_string_tuple(
                self.excluded_endpoint_ids,
                "PredictionFirstOutcomeContract.excluded_endpoint_ids",
                allow_empty=True,
            ),
        )
        for field_name in (
            "forbid_fit",
            "forbid_model_selection",
            "require_prediction_commit",
        ):
            _set_frozen(
                self,
                field_name,
                _as_bool(
                    getattr(self, field_name),
                    f"PredictionFirstOutcomeContract.{field_name}",
                ),
            )
        if self.schema_version != "masld-bench-prediction-first-outcome-v1":
            raise ContractError("unsupported prediction-first outcome schema_version")
        expected = {
            "scope": "non_champion_stress_only",
            "design_visibility": "development_visible",
            "outcome_visibility": "withheld_until_prediction_commit",
            "prediction_scope": "genome_wide",
        }
        for field_name, expected_value in expected.items():
            if getattr(self, field_name) != expected_value:
                raise ContractError(
                    f"PredictionFirstOutcomeContract.{field_name} must be "
                    f"{expected_value!r}"
                )
        if not (
            self.forbid_fit
            and self.forbid_model_selection
            and self.require_prediction_commit
        ):
            raise ContractError(
                "prediction-first outcomes must forbid fitting/model selection and "
                "require a frozen prediction commit"
            )
        if (
            len(set(self.excluded_endpoint_ids)) != len(self.excluded_endpoint_ids)
            or tuple(sorted(self.excluded_endpoint_ids))
            != self.excluded_endpoint_ids
        ):
            raise ContractError(
                "prediction-first excluded_endpoint_ids must be unique and sorted"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PredictionFirstOutcomeContract":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS,
            label="PredictionFirstOutcomeContract",
        )
        return cls(**raw)


@dataclass(frozen=True, slots=True)
class DatasetManifest(StrictContract):
    schema_version: str
    dataset_id: str
    title: str
    role: DatasetRole
    status: DatasetStatus
    species: str
    accession: tuple[str, ...]
    access_tier: AccessState
    automatic_download: bool
    outcome_url: str | None
    redistribution_class: str
    license_status: str
    biological_unit: str
    expected_biological_units: int | str | None
    native_genome_build: str
    native_annotation_release: str
    modalities: tuple[str, ...]
    pairing_levels: tuple[PairingState, ...]
    modality_status_default: MissingState
    exposure_status: ExposureState
    admission_blocking: bool
    blockers: tuple[str, ...]
    notes: str
    provenance: DatasetProvenance
    split: DatasetSplit
    activation: DatasetActivationContract = field(
        default_factory=DatasetActivationContract
    )
    prediction_first_outcome: PredictionFirstOutcomeContract | None = None

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "dataset_id",
            "title",
            "role",
            "status",
            "species",
            "accession",
            "access_tier",
            "automatic_download",
            "outcome_url",
            "redistribution_class",
            "license_status",
            "biological_unit",
            "expected_biological_units",
            "native_genome_build",
            "native_annotation_release",
            "modalities",
            "pairing_levels",
            "modality_status_default",
            "exposure_status",
            "admission_blocking",
            "blockers",
            "notes",
            "provenance",
            "split",
            "activation",
            "prediction_first_outcome",
        }
    )
    _OPTIONAL: ClassVar[frozenset[str]] = frozenset(
        {
            "outcome_url",
            "expected_biological_units",
            "blockers",
            "notes",
            "activation",
            "prediction_first_outcome",
        }
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "dataset_id",
            "title",
            "species",
            "redistribution_class",
            "license_status",
            "biological_unit",
            "native_genome_build",
            "native_annotation_release",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"DatasetManifest.{field_name}"),
            )
        _set_frozen(self, "role", _as_enum(self.role, DatasetRole, "DatasetManifest.role"))
        _set_frozen(
            self, "status", _as_enum(self.status, DatasetStatus, "DatasetManifest.status")
        )
        _set_frozen(
            self,
            "accession",
            _as_string_tuple(
                self.accession, "DatasetManifest.accession", allow_empty=False
            ),
        )
        _set_frozen(
            self,
            "access_tier",
            _as_enum(self.access_tier, AccessState, "DatasetManifest.access_tier"),
        )
        _set_frozen(
            self,
            "automatic_download",
            _as_bool(self.automatic_download, "DatasetManifest.automatic_download"),
        )
        _set_frozen(
            self,
            "outcome_url",
            _optional_string(self.outcome_url, "DatasetManifest.outcome_url"),
        )
        _set_frozen(
            self,
            "expected_biological_units",
            _optional_int_or_unresolved(
                self.expected_biological_units,
                "DatasetManifest.expected_biological_units",
                minimum=1,
            ),
        )
        _set_frozen(
            self,
            "modalities",
            _as_string_tuple(
                self.modalities, "DatasetManifest.modalities", allow_empty=False
            ),
        )
        _set_frozen(
            self,
            "pairing_levels",
            _as_enum_tuple(
                self.pairing_levels,
                PairingState,
                "DatasetManifest.pairing_levels",
                allow_empty=False,
            ),
        )
        _set_frozen(
            self,
            "modality_status_default",
            _as_enum(
                self.modality_status_default,
                MissingState,
                "DatasetManifest.modality_status_default",
            ),
        )
        _set_frozen(
            self,
            "exposure_status",
            _as_enum(
                self.exposure_status,
                ExposureState,
                "DatasetManifest.exposure_status",
            ),
        )
        _set_frozen(
            self,
            "admission_blocking",
            _as_bool(self.admission_blocking, "DatasetManifest.admission_blocking"),
        )
        _set_frozen(
            self, "blockers", _as_string_tuple(self.blockers, "DatasetManifest.blockers")
        )
        is_withheld_role = self.role is DatasetRole.WITHHELD_SEALED
        is_withheld_status = self.status is DatasetStatus.WITHHELD_SEALED
        if is_withheld_role != is_withheld_status:
            raise ContractError(
                "withheld datasets must use both role=withheld_sealed and "
                "status=withheld_sealed"
            )
        if is_withheld_role:
            if self.automatic_download:
                raise ContractError("withheld datasets may not enable automatic download")
            if self.outcome_url != "WITHHELD_SEALED":
                raise ContractError(
                    "withheld dataset outcome_url must be the WITHHELD_SEALED sentinel"
                )
            if self.split.label_visibility != "withheld_sealed":
                raise ContractError("withheld dataset labels must remain withheld_sealed")
            if self.modality_status_default is not MissingState.WITHHELD_SEALED:
                raise ContractError(
                    "withheld datasets must default modalities to withheld_sealed"
                )
            if not self.admission_blocking or not self.blockers:
                raise ContractError(
                    "withheld datasets remain admission-blocking until an evaluator activates them"
                )
        if self.status in {DatasetStatus.BLOCKED, DatasetStatus.DEFERRED} and (
            not self.admission_blocking or not self.blockers
        ):
            raise ContractError(
                "blocked or deferred datasets require admission_blocking=true and blockers"
            )
        _set_frozen(
            self,
            "notes",
            _as_string(self.notes, "DatasetManifest.notes", allow_empty=True),
        )
        if not isinstance(self.provenance, DatasetProvenance):
            raise ContractError("DatasetManifest.provenance must be DatasetProvenance")
        if not isinstance(self.split, DatasetSplit):
            raise ContractError("DatasetManifest.split must be DatasetSplit")
        if not isinstance(self.activation, DatasetActivationContract):
            raise ContractError(
                "DatasetManifest.activation must be DatasetActivationContract"
            )
        if self.prediction_first_outcome is not None:
            if not isinstance(
                self.prediction_first_outcome, PredictionFirstOutcomeContract
            ):
                raise ContractError(
                    "DatasetManifest.prediction_first_outcome must be a "
                    "PredictionFirstOutcomeContract or null"
                )
            if self.role is not DatasetRole.MECHANISTIC_DIAGNOSTIC:
                raise ContractError(
                    "prediction-first outcomes are restricted to mechanistic diagnostics"
                )
            if self.split.label_visibility != "prediction_first_outcomes_withheld":
                raise ContractError(
                    "prediction-first dataset split must keep outcomes withheld"
                )
            if self.outcome_url != "PREDICTION_FIRST_WITHHELD":
                raise ContractError(
                    "prediction-first dataset outcome_url must use the withheld sentinel"
                )
            expected_universe_role = f"prediction_universe:{self.dataset_id}"
            if (
                self.prediction_first_outcome.prediction_universe_role
                != expected_universe_role
            ):
                raise ContractError(
                    "prediction-first universe role must bind the DatasetManifest.dataset_id"
                )
        expected_evidence_roles = {
            f"dataset_authority:{self.dataset_id}:{authority}"
            for authority in DatasetActivationContract._AUTHORITY_FIELDS
        }
        observed_evidence_roles = {
            artifact.role for artifact in self.activation.evidence_artifacts
        }
        if observed_evidence_roles and observed_evidence_roles != expected_evidence_roles:
            raise ContractError(
                "dataset activation evidence roles must bind the DatasetManifest.dataset_id"
            )
        if not self.admission_blocking and not self.activation.ready:
            raise ContractError(
                "an admission-ready dataset must bind a complete DatasetActivationContract"
            )
        if (
            self.exposure_status is ExposureState.UNKNOWN
            and not self.admission_blocking
        ):
            raise ContractError("unknown dataset exposure must be admission-blocking")
        if "UNRESOLVED" in self.native_genome_build.upper() and not self.admission_blocking:
            raise ContractError("unresolved native genome build must be admission-blocking")
        if (
            "UNRESOLVED" in self.native_annotation_release.upper()
            and not self.admission_blocking
        ):
            raise ContractError(
                "unresolved native annotation release must be admission-blocking"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DatasetManifest":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS - cls._OPTIONAL,
            label="DatasetManifest",
        )
        return cls(
            schema_version=raw["schema_version"],
            dataset_id=raw["dataset_id"],
            title=raw["title"],
            role=raw["role"],
            status=raw["status"],
            species=raw["species"],
            accession=raw["accession"],
            access_tier=raw["access_tier"],
            automatic_download=raw["automatic_download"],
            outcome_url=raw.get("outcome_url"),
            redistribution_class=raw["redistribution_class"],
            license_status=raw["license_status"],
            biological_unit=raw["biological_unit"],
            expected_biological_units=raw.get("expected_biological_units"),
            native_genome_build=raw["native_genome_build"],
            native_annotation_release=raw["native_annotation_release"],
            modalities=raw["modalities"],
            pairing_levels=raw["pairing_levels"],
            modality_status_default=raw["modality_status_default"],
            exposure_status=raw["exposure_status"],
            admission_blocking=raw["admission_blocking"],
            blockers=raw.get("blockers", ()),
            notes=raw.get("notes", ""),
            provenance=DatasetProvenance.from_dict(raw["provenance"]),
            split=DatasetSplit.from_dict(raw["split"]),
            activation=(
                DatasetActivationContract.from_dict(raw["activation"])
                if "activation" in raw
                else DatasetActivationContract()
            ),
            prediction_first_outcome=(
                PredictionFirstOutcomeContract.from_dict(
                    raw["prediction_first_outcome"]
                )
                if raw.get("prediction_first_outcome") is not None
                else None
            ),
        )


@dataclass(frozen=True, slots=True)
class ModelExecutionContract(StrictContract):
    """Exact repository, preprocessing, terms, runtime, and adapter binding."""

    schema_version: str = "masld-bench-model-execution-v1"
    upstream_repository: str = "UNRESOLVED"
    upstream_commit: str = "UNRESOLVED"
    preprocessing_sha256: str = "UNRESOLVED"
    vocabulary_or_build: str = "UNRESOLVED"
    genome_build: str = "UNRESOLVED"
    code_license: str = "UNRESOLVED"
    weights_license: str = "UNRESOLVED"
    derivative_weights_license: str = "UNRESOLVED"
    training_cutoff: str = "UNRESOLVED"
    declared_corpora_sha256: str = "UNRESOLVED"
    cellxgene_uuids_sha256: str = "UNRESOLVED"
    exposure_audit_sha256: str = "UNRESOLVED"
    evidence_artifacts: tuple[ArtifactRef, ...] = ()
    supported_actions: tuple[str, ...] = ()
    executable_tasks: tuple[str, ...] = ()
    supported_adaptation_regimes: tuple[str, ...] = ()
    runtime_id: str = "UNRESOLVED"
    adapter_command: tuple[str, ...] = ()
    environment_artifact: ArtifactRef | None = None
    checkpoint_artifact: ArtifactRef | None = None
    capture_r_session: bool = False
    ready: bool = False

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "upstream_repository",
            "upstream_commit",
            "preprocessing_sha256",
            "vocabulary_or_build",
            "genome_build",
            "code_license",
            "weights_license",
            "derivative_weights_license",
            "training_cutoff",
            "declared_corpora_sha256",
            "cellxgene_uuids_sha256",
            "exposure_audit_sha256",
            "evidence_artifacts",
            "supported_actions",
            "executable_tasks",
            "supported_adaptation_regimes",
            "runtime_id",
            "adapter_command",
            "environment_artifact",
            "checkpoint_artifact",
            "capture_r_session",
            "ready",
        }
    )
    _OPTIONAL: ClassVar[frozenset[str]] = frozenset(
        {
            "evidence_artifacts",
            "executable_tasks",
            "supported_adaptation_regimes",
            "environment_artifact",
            "checkpoint_artifact",
        }
    )
    _ACTIONS: ClassVar[frozenset[str]] = frozenset(
        {"probe", "prepare", "fit", "predict", "export"}
    )
    _AUTHORITY_FIELDS: ClassVar[tuple[str, ...]] = (
        "preprocessing",
        "code_license",
        "weights_license",
        "derivative_weights_license",
        "declared_corpora",
        "cellxgene_uuids",
        "exposure_audit",
    )
    _KNOWN_CODE_LICENSES: ClassVar[frozenset[str]] = frozenset(
        {
            "MIT",
            "BSD-2-Clause",
            "BSD-3-Clause",
            "Apache-2.0",
            "GPL-2.0-only",
            "GPL-2.0-or-later",
            "GPL-3.0-only",
            "GPL-3.0-or-later",
            "LGPL-3.0-only",
            "LGPL-3.0-or-later",
            "project_owned",
            "research_only",
            "noncommercial",
            "restricted",
            "source_terms_apply",
        }
    )
    _KNOWN_WEIGHT_LICENSES: ClassVar[frozenset[str]] = frozenset(
        {
            "MIT",
            "BSD-2-Clause",
            "BSD-3-Clause",
            "Apache-2.0",
            "CC-BY-4.0",
            "CC-BY-SA-4.0",
            "CC-BY-NC-4.0",
            "CC-BY-NC-SA-4.0",
            "CC-BY-NC-ND-4.0",
            "CC0-1.0",
            "GPL-2.0-only",
            "GPL-2.0-or-later",
            "GPL-3.0-only",
            "GPL-3.0-or-later",
            "LGPL-3.0-only",
            "LGPL-3.0-or-later",
            "project_owned",
            "not_applicable_classical",
            "research_only",
            "noncommercial",
            "restricted",
            "redistribution_prohibited",
            "source_terms_apply",
        }
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "upstream_repository",
            "upstream_commit",
            "vocabulary_or_build",
            "genome_build",
            "code_license",
            "weights_license",
            "derivative_weights_license",
            "training_cutoff",
            "runtime_id",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(
                    getattr(self, field_name),
                    f"ModelExecutionContract.{field_name}",
                ),
            )
        if self.schema_version != "masld-bench-model-execution-v1":
            raise ContractError(
                "ModelExecutionContract.schema_version must be "
                "masld-bench-model-execution-v1"
            )
        for field_name in (
            "preprocessing_sha256",
            "declared_corpora_sha256",
            "cellxgene_uuids_sha256",
            "exposure_audit_sha256",
        ):
            _set_frozen(
                self,
                field_name,
                _as_sha256_or_unresolved(
                    getattr(self, field_name),
                    f"ModelExecutionContract.{field_name}",
                ),
            )
        actions = _as_string_tuple(
            self.supported_actions,
            "ModelExecutionContract.supported_actions",
        )
        unknown_actions = sorted(set(actions).difference(self._ACTIONS))
        if unknown_actions:
            raise ContractError(
                "ModelExecutionContract has unsupported action(s): "
                + ", ".join(unknown_actions)
            )
        _set_frozen(self, "supported_actions", actions)
        executable_tasks = _as_string_tuple(
            self.executable_tasks,
            "ModelExecutionContract.executable_tasks",
        )
        if executable_tasks != tuple(sorted(set(executable_tasks))):
            raise ContractError(
                "ModelExecutionContract.executable_tasks must be sorted and unique"
            )
        _set_frozen(self, "executable_tasks", executable_tasks)
        adaptation_regimes = _as_string_tuple(
            self.supported_adaptation_regimes,
            "ModelExecutionContract.supported_adaptation_regimes",
        )
        if any(not _IDENTIFIER_PATTERN.fullmatch(item) for item in adaptation_regimes):
            raise ContractError(
                "ModelExecutionContract.supported_adaptation_regimes must contain "
                "lowercase identifiers"
            )
        _set_frozen(
            self,
            "supported_adaptation_regimes",
            adaptation_regimes,
        )
        command = _as_string_tuple(
            self.adapter_command,
            "ModelExecutionContract.adapter_command",
        )
        _set_frozen(self, "adapter_command", command)
        if self.environment_artifact is not None and not isinstance(
            self.environment_artifact, ArtifactRef
        ):
            raise ContractError(
                "ModelExecutionContract.environment_artifact must be ArtifactRef or null"
            )
        if self.checkpoint_artifact is not None and not isinstance(
            self.checkpoint_artifact, ArtifactRef
        ):
            raise ContractError(
                "ModelExecutionContract.checkpoint_artifact must be ArtifactRef or null"
            )
        if not isinstance(self.evidence_artifacts, (tuple, list)):
            raise ContractError("ModelExecutionContract.evidence_artifacts must be an array")
        evidence_artifacts = tuple(self.evidence_artifacts)
        if any(not isinstance(item, ArtifactRef) for item in evidence_artifacts):
            raise ContractError(
                "ModelExecutionContract.evidence_artifacts entries must be ArtifactRef"
            )
        evidence_paths = [item.path for item in evidence_artifacts]
        evidence_roles = [item.role for item in evidence_artifacts]
        if len(set(evidence_paths)) != len(evidence_paths):
            raise ContractError(
                "ModelExecutionContract.evidence_artifacts paths must be unique"
            )
        if any(role is None for role in evidence_roles) or len(set(evidence_roles)) != len(
            evidence_roles
        ):
            raise ContractError(
                "ModelExecutionContract.evidence_artifacts require unique roles"
            )
        suffix_to_artifact: dict[str, ArtifactRef] = {}
        for artifact in evidence_artifacts:
            assert artifact.role is not None
            parts = artifact.role.split(":")
            if (
                len(parts) != 3
                or parts[0] != "model_authority"
                or not parts[1]
                or parts[2] not in self._AUTHORITY_FIELDS
            ):
                raise ContractError(
                    "ModelExecutionContract evidence role must be "
                    "model_authority:<model_id>:<authority>"
                )
            suffix_to_artifact[parts[2]] = artifact
        if evidence_artifacts and set(suffix_to_artifact) != set(self._AUTHORITY_FIELDS):
            raise ContractError(
                "ModelExecutionContract.evidence_artifacts must cover every and only "
                "preprocessing, code_license, weights_license, "
                "derivative_weights_license, declared_corpora, cellxgene_uuids, "
                "and exposure_audit"
            )
        for authority in (
            "preprocessing",
            "declared_corpora",
            "cellxgene_uuids",
            "exposure_audit",
        ):
            artifact = suffix_to_artifact.get(authority)
            if artifact is not None and artifact.sha256 != getattr(
                self, f"{authority}_sha256"
            ):
                raise ContractError(
                    "ModelExecutionContract authority artifact SHA-256 differs from "
                    f"{authority}_sha256"
                )
        _set_frozen(self, "evidence_artifacts", evidence_artifacts)
        _set_frozen(
            self,
            "capture_r_session",
            _as_bool(
                self.capture_r_session,
                "ModelExecutionContract.capture_r_session",
            ),
        )
        _set_frozen(
            self,
            "ready",
            _as_bool(self.ready, "ModelExecutionContract.ready"),
        )
        if self.ready:
            if self.code_license not in self._KNOWN_CODE_LICENSES:
                raise ContractError(
                    "ready model execution has an unrecognized code_license"
                )
            for field_name in ("weights_license", "derivative_weights_license"):
                if getattr(self, field_name) not in self._KNOWN_WEIGHT_LICENSES:
                    raise ContractError(
                        f"ready model execution has an unrecognized {field_name}"
                    )
        unresolved = (
            any(
                _contains_unresolved_marker(getattr(self, field_name))
                for field_name in (
                    "upstream_repository",
                    "upstream_commit",
                    "preprocessing_sha256",
                    "vocabulary_or_build",
                    "genome_build",
                    "code_license",
                    "weights_license",
                    "derivative_weights_license",
                    "training_cutoff",
                    "declared_corpora_sha256",
                    "cellxgene_uuids_sha256",
                    "exposure_audit_sha256",
                    "runtime_id",
                )
            )
            or not actions
            or not adaptation_regimes
            or not command
            or self.environment_artifact is None
            or len(evidence_artifacts) != len(self._AUTHORITY_FIELDS)
            or any(_contains_unresolved_marker(part) for part in command)
        )
        if self.ready and unresolved:
            raise ContractError(
                "ready model execution requires exact repository, preprocessing, terms, "
                "exposure, runtime, supported actions, adaptation regimes, and adapter "
                "command"
            )
        if not self.ready and not unresolved:
            raise ContractError(
                "complete model execution evidence must set ready=true"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ModelExecutionContract":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS - cls._OPTIONAL,
            label="ModelExecutionContract",
        )
        environment = raw.get("environment_artifact")
        checkpoint = raw.get("checkpoint_artifact")
        evidence = raw.get("evidence_artifacts", [])
        return cls(
            schema_version=raw["schema_version"],
            upstream_repository=raw["upstream_repository"],
            upstream_commit=raw["upstream_commit"],
            preprocessing_sha256=raw["preprocessing_sha256"],
            vocabulary_or_build=raw["vocabulary_or_build"],
            genome_build=raw["genome_build"],
            code_license=raw["code_license"],
            weights_license=raw["weights_license"],
            derivative_weights_license=raw["derivative_weights_license"],
            training_cutoff=raw["training_cutoff"],
            declared_corpora_sha256=raw["declared_corpora_sha256"],
            cellxgene_uuids_sha256=raw["cellxgene_uuids_sha256"],
            exposure_audit_sha256=raw["exposure_audit_sha256"],
            evidence_artifacts=tuple(
                ArtifactRef.from_dict(item) for item in evidence
            ),
            supported_actions=raw["supported_actions"],
            executable_tasks=raw.get("executable_tasks", ()),
            supported_adaptation_regimes=raw.get(
                "supported_adaptation_regimes", ()
            ),
            runtime_id=raw["runtime_id"],
            adapter_command=raw["adapter_command"],
            environment_artifact=(
                ArtifactRef.from_dict(environment) if environment is not None else None
            ),
            checkpoint_artifact=(
                ArtifactRef.from_dict(checkpoint) if checkpoint is not None else None
            ),
            capture_r_session=raw["capture_r_session"],
            ready=raw["ready"],
        )


@dataclass(frozen=True, slots=True)
class ModelManifest(StrictContract):
    schema_version: str
    family_id: str
    modality: tuple[str, ...]
    family_status: ReleaseClass
    model_id: str
    display_name: str
    implementation_type: ImplementationType
    checkpoint_revision: str
    checkpoint_sha256: str
    license_status: str
    exposure_status: ExposureState
    status: ReleaseClass
    admission_blocking: bool
    blockers: tuple[str, ...]
    supported_tasks: tuple[str, ...]
    execution: ModelExecutionContract = field(default_factory=ModelExecutionContract)

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "family_id",
            "modality",
            "family_status",
            "model_id",
            "display_name",
            "implementation_type",
            "checkpoint_revision",
            "checkpoint_sha256",
            "license_status",
            "exposure_status",
            "status",
            "admission_blocking",
            "blockers",
            "supported_tasks",
            "execution",
        }
    )
    _OPTIONAL: ClassVar[frozenset[str]] = frozenset({"execution"})
    _FAMILY_FIELDS: ClassVar[frozenset[str]] = frozenset(
        {"schema_version", "family_id", "modality", "status", "models"}
    )
    _ENTRY_FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "model_id",
            "display_name",
            "implementation_type",
            "checkpoint_revision",
            "checkpoint_sha256",
            "license_status",
            "exposure_status",
            "status",
            "admission_blocking",
            "blockers",
            "supported_tasks",
            "execution",
        }
    )
    _ENTRY_OPTIONAL: ClassVar[frozenset[str]] = frozenset({"execution"})

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "family_id",
            "model_id",
            "display_name",
            "checkpoint_revision",
            "license_status",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"ModelManifest.{field_name}"),
            )
        _set_frozen(
            self,
            "modality",
            _as_string_tuple(
                self.modality, "ModelManifest.modality", allow_empty=False
            ),
        )
        _set_frozen(
            self,
            "family_status",
            _as_enum(
                self.family_status, ReleaseClass, "ModelManifest.family_status"
            ),
        )
        _set_frozen(
            self,
            "implementation_type",
            _as_enum(
                self.implementation_type,
                ImplementationType,
                "ModelManifest.implementation_type",
            ),
        )
        _set_frozen(
            self,
            "checkpoint_sha256",
            _as_sha256_or_unresolved(
                self.checkpoint_sha256, "ModelManifest.checkpoint_sha256"
            ),
        )
        _set_frozen(
            self,
            "exposure_status",
            _as_enum(
                self.exposure_status, ExposureState, "ModelManifest.exposure_status"
            ),
        )
        _set_frozen(
            self, "status", _as_enum(self.status, ReleaseClass, "ModelManifest.status")
        )
        _set_frozen(
            self,
            "admission_blocking",
            _as_bool(self.admission_blocking, "ModelManifest.admission_blocking"),
        )
        _set_frozen(
            self, "blockers", _as_string_tuple(self.blockers, "ModelManifest.blockers")
        )
        _set_frozen(
            self,
            "supported_tasks",
            _as_string_tuple(
                self.supported_tasks,
                "ModelManifest.supported_tasks",
                allow_empty=True,
            ),
        )
        if not self.supported_tasks and self.status is not ReleaseClass.DEFERRED:
            raise ContractError(
                "only deferred models may have no supported benchmark tasks"
            )
        unresolved = (
            any(
                "UNRESOLVED" in value.upper()
                for value in (
                    self.checkpoint_revision,
                    self.checkpoint_sha256,
                    self.license_status,
                )
            )
            or self.exposure_status is ExposureState.UNKNOWN
        )
        if unresolved and not self.admission_blocking:
            raise ContractError("UNRESOLVED model fields must be admission-blocking")
        if (
            self.checkpoint_sha256 == "NOT_APPLICABLE"
            and self.implementation_type is ImplementationType.PRETRAINED_ADAPTER
        ):
            raise ContractError(
                "pretrained adapters require an exact checkpoint SHA-256 or an admission-blocking UNRESOLVED value"
            )
        if not isinstance(self.execution, ModelExecutionContract):
            raise ContractError(
                "ModelManifest.execution must be ModelExecutionContract"
            )
        if not self.admission_blocking and not self.execution.ready:
            raise ContractError(
                "an admission-ready model must bind a complete ModelExecutionContract"
            )
        if self.execution.ready:
            if self.execution.executable_tasks and not set(
                self.execution.executable_tasks
            ).issubset(self.supported_tasks):
                raise ContractError(
                    "model execution executable_tasks must be a subset of supported_tasks"
                )
            expected_evidence_roles = {
                f"model_authority:{self.model_id}:{authority}"
                for authority in ModelExecutionContract._AUTHORITY_FIELDS
            }
            observed_evidence_roles = {
                artifact.role for artifact in self.execution.evidence_artifacts
            }
            if observed_evidence_roles != expected_evidence_roles:
                raise ContractError(
                    "model execution evidence roles must bind the ModelManifest.model_id"
                )
            checkpoint = self.execution.checkpoint_artifact
            if self.implementation_type is ImplementationType.PRETRAINED_ADAPTER:
                if checkpoint is None:
                    raise ContractError(
                        "ready pretrained adapters require a checkpoint ArtifactRef"
                    )
                if checkpoint.sha256 != self.checkpoint_sha256:
                    raise ContractError(
                        "model checkpoint SHA-256 must match its checkpoint ArtifactRef"
                    )
            elif checkpoint is not None and self.checkpoint_sha256 == "NOT_APPLICABLE":
                raise ContractError(
                    "models with checkpoint_sha256=NOT_APPLICABLE may not bind a checkpoint artifact"
                )

    @property
    def family(self) -> str:
        return self.family_id

    @property
    def title(self) -> str:
        return self.display_name

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ModelManifest":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS - cls._OPTIONAL,
            label="ModelManifest",
        )
        normalized = dict(raw)
        normalized["execution"] = (
            ModelExecutionContract.from_dict(raw["execution"])
            if "execution" in raw
            else ModelExecutionContract()
        )
        return cls(**normalized)

    @classmethod
    def from_family_entry(
        cls,
        family: Mapping[str, Any],
        model: Mapping[str, Any],
    ) -> "ModelManifest":
        family_raw = _strict_mapping(
            family,
            allowed=cls._FAMILY_FIELDS,
            required=cls._FAMILY_FIELDS,
            label="ModelFamilyManifest",
        )
        model_raw = _strict_mapping(
            model,
            allowed=cls._ENTRY_FIELDS,
            required=cls._ENTRY_FIELDS - cls._ENTRY_OPTIONAL,
            label=f"ModelFamilyManifest.models[{model.get('model_id', '?') if isinstance(model, Mapping) else '?'}]",
        )
        normalized_model = dict(model_raw)
        normalized_model["execution"] = (
            ModelExecutionContract.from_dict(model_raw["execution"])
            if "execution" in model_raw
            else ModelExecutionContract()
        )
        return cls(
            schema_version=family_raw["schema_version"],
            family_id=family_raw["family_id"],
            modality=family_raw["modality"],
            family_status=family_raw["status"],
            **normalized_model,
        )


@dataclass(frozen=True, slots=True)
class LocusSplitControls(StrictContract):
    """Executable anti-leakage controls for locus-level outer splits."""

    schema_version: str
    block_merge_rule: str
    merge_distance_bp: int
    ld_r2_threshold: float
    ld_ancestry_rule: str
    sequence_identity_grouping: tuple[str, ...]
    boundary_buffer_policy: str

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "block_merge_rule",
            "merge_distance_bp",
            "ld_r2_threshold",
            "ld_ancestry_rule",
            "sequence_identity_grouping",
            "boundary_buffer_policy",
        }
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "block_merge_rule",
            "ld_ancestry_rule",
            "boundary_buffer_policy",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(
                    getattr(self, field_name),
                    f"LocusSplitControls.{field_name}",
                ),
            )
        _set_frozen(
            self,
            "merge_distance_bp",
            _as_int(
                self.merge_distance_bp,
                "LocusSplitControls.merge_distance_bp",
                minimum=1,
            ),
        )
        if (
            isinstance(self.ld_r2_threshold, bool)
            or not isinstance(self.ld_r2_threshold, (int, float))
        ):
            raise ContractError(
                "LocusSplitControls.ld_r2_threshold must be numeric"
            )
        _set_frozen(self, "ld_r2_threshold", float(self.ld_r2_threshold))
        _set_frozen(
            self,
            "sequence_identity_grouping",
            _as_string_tuple(
                self.sequence_identity_grouping,
                "LocusSplitControls.sequence_identity_grouping",
                allow_empty=False,
            ),
        )
        expected = {
            "schema_version": "masld-bench-locus-split-controls-v1",
            "block_merge_rule": "within_distance_or_ld_threshold",
            "merge_distance_bp": 1_000_000,
            "ld_r2_threshold": 0.8,
            "ld_ancestry_rule": "any_represented_ancestry",
            "sequence_identity_grouping": (
                "exact_forward_identity",
                "exact_reverse_complement_identity",
            ),
            "boundary_buffer_policy": (
                "exclude_each_outer_boundary_by_largest_admitted_model_receptive_field"
            ),
        }
        for field_name, expected_value in expected.items():
            if getattr(self, field_name) != expected_value:
                raise ContractError(
                    f"LocusSplitControls.{field_name} must be {expected_value!r}"
                )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "LocusSplitControls":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS,
            label="LocusSplitControls",
        )
        return cls(**raw)


@dataclass(frozen=True, slots=True)
class SplitSpec(StrictContract):
    """Strict grouped-split identity bound into every benchmark run."""

    schema_version: str
    split_id: str
    entity: str
    group_keys: tuple[str, ...]
    roles: tuple[str, ...]
    seed: int
    locus_grouping: str
    label_policy: str
    exposure_policy: str
    admission_gates: tuple[str, ...]
    locus_controls: LocusSplitControls | None = None

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "split_id",
            "entity",
            "group_keys",
            "roles",
            "seed",
            "locus_grouping",
            "label_policy",
            "exposure_policy",
            "admission_gates",
            "locus_controls",
        }
    )
    _OPTIONAL: ClassVar[frozenset[str]] = frozenset({"locus_controls"})

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "split_id",
            "entity",
            "locus_grouping",
            "label_policy",
            "exposure_policy",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"SplitSpec.{field_name}"),
            )
        for field_name in ("group_keys", "roles", "admission_gates"):
            _set_frozen(
                self,
                field_name,
                _as_string_tuple(
                    getattr(self, field_name),
                    f"SplitSpec.{field_name}",
                    allow_empty=False,
                ),
            )
        _set_frozen(self, "seed", _as_int(self.seed, "SplitSpec.seed"))
        if self.schema_version != "masld-bench-split-v1":
            raise ContractError("SplitSpec.schema_version is unsupported")
        if self.roles != ("train", "development", "withheld_sealed"):
            raise ContractError(
                "SplitSpec.roles must be train, development, withheld_sealed in order"
            )
        if self.seed != 20260821:
            raise ContractError("SplitSpec.seed must be the frozen seed 20260821")
        if self.entity in {"cell", "spot", "section", "array", "run", "guide"}:
            raise ContractError(
                "SplitSpec.entity must be a biological donor/participant or locus unit"
            )
        if self.split_id == "locus_outer":
            if not isinstance(self.locus_controls, LocusSplitControls):
                raise ContractError(
                    "locus_outer requires strict LocusSplitControls"
                )
            if "sequence_identity_group" not in self.group_keys:
                raise ContractError(
                    "locus_outer must group sequence-identity equivalents"
                )
        elif self.locus_controls is not None:
            raise ContractError(
                "locus_controls are permitted only for locus_outer"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SplitSpec":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS - cls._OPTIONAL,
            label="SplitSpec",
        )
        return cls(
            schema_version=raw["schema_version"],
            split_id=raw["split_id"],
            entity=raw["entity"],
            group_keys=raw["group_keys"],
            roles=raw["roles"],
            seed=raw["seed"],
            locus_grouping=raw["locus_grouping"],
            label_policy=raw["label_policy"],
            exposure_policy=raw["exposure_policy"],
            admission_gates=raw["admission_gates"],
            locus_controls=(
                LocusSplitControls.from_dict(raw["locus_controls"])
                if "locus_controls" in raw
                else None
            ),
        )


@dataclass(frozen=True, slots=True)
class TaskSpec(StrictContract):
    schema_version: str
    task_id: str
    status: TaskStatus
    unit_of_inference: str
    input_modalities: tuple[str, ...]
    endpoint: str
    metrics: tuple[str, ...]
    datasets_train: tuple[str, ...]
    datasets_development: tuple[str, ...]
    datasets_sealed: tuple[str, ...]
    split_id: str
    required_pairing_levels: tuple[PairingState, ...]
    missingness_policy: str
    admission_gates: tuple[str, ...]
    baseline_model_ids: tuple[str, ...]
    uncertainty_method: str
    resampling_units: tuple[str, ...]
    bootstrap_replicates: int
    multiplicity_family: str
    multiplicity_method: str
    promotion_gate_id: str
    promotion_gate_config_sha256: str
    primary_evaluator_id: str
    evaluator_parameters: Mapping[str, Any]
    claim_gate: str

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "task_id",
            "status",
            "unit_of_inference",
            "input_modalities",
            "endpoint",
            "metrics",
            "datasets_train",
            "datasets_development",
            "datasets_sealed",
            "split_id",
            "required_pairing_levels",
            "missingness_policy",
            "admission_gates",
            "baseline_model_ids",
            "uncertainty_method",
            "resampling_units",
            "bootstrap_replicates",
            "multiplicity_family",
            "multiplicity_method",
            "promotion_gate_id",
            "promotion_gate_config_sha256",
            "primary_evaluator_id",
            "evaluator_parameters",
            "claim_gate",
        }
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "task_id",
            "unit_of_inference",
            "endpoint",
            "split_id",
            "missingness_policy",
            "uncertainty_method",
            "multiplicity_family",
            "multiplicity_method",
            "promotion_gate_id",
            "primary_evaluator_id",
            "claim_gate",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"TaskSpec.{field_name}"),
            )
        _set_frozen(self, "status", _as_enum(self.status, TaskStatus, "TaskSpec.status"))
        for field_name, allow_empty in (
            ("input_modalities", False),
            ("metrics", False),
            ("datasets_train", True),
            ("datasets_development", True),
            ("datasets_sealed", True),
            ("admission_gates", True),
            ("baseline_model_ids", False),
            ("resampling_units", False),
        ):
            _set_frozen(
                self,
                field_name,
                _as_string_tuple(
                    getattr(self, field_name),
                    f"TaskSpec.{field_name}",
                    allow_empty=allow_empty,
                ),
            )
        _set_frozen(
            self,
            "bootstrap_replicates",
            _as_int(
                self.bootstrap_replicates,
                "TaskSpec.bootstrap_replicates",
                minimum=1,
            ),
        )
        _set_frozen(
            self,
            "promotion_gate_config_sha256",
            _as_sha256(
                self.promotion_gate_config_sha256,
                "TaskSpec.promotion_gate_config_sha256",
            ),
        )
        _set_frozen(
            self,
            "evaluator_parameters",
            _as_metadata(
                self.evaluator_parameters,
                "TaskSpec.evaluator_parameters",
            ),
        )
        _set_frozen(
            self,
            "required_pairing_levels",
            _as_enum_tuple(
                self.required_pairing_levels,
                PairingState,
                "TaskSpec.required_pairing_levels",
                allow_empty=True,
            ),
        )
        groups = (
            set(self.datasets_train),
            set(self.datasets_development),
            set(self.datasets_sealed),
        )
        if groups[0] & groups[1] or groups[0] & groups[2] or groups[1] & groups[2]:
            raise ContractError("TaskSpec dataset partitions must be disjoint")

    @property
    def dataset_ids(self) -> tuple[str, ...]:
        return (
            *self.datasets_train,
            *self.datasets_development,
            *self.datasets_sealed,
        )

    @property
    def primary_metric(self) -> str:
        return self.metrics[0]

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "TaskSpec":
        raw = _strict_mapping(
            value, allowed=cls._FIELDS, required=cls._FIELDS, label="TaskSpec"
        )
        return cls(**raw)


@dataclass(frozen=True, slots=True)
class RunSpec(StrictContract):
    schema_version: str
    campaign_id: str
    task_id: str
    model_id: str
    dataset_ids: tuple[str, ...]
    split_id: str
    seed: int
    stage: str = "benchmark"
    adaptation_regime: str = "common_lane"
    fold: int = 0
    action: tuple[str, ...] = ()
    adapter_command: tuple[str, ...] = ()
    immutable_inputs: Mapping[str, str] = field(default_factory=dict)
    inputs: tuple[ArtifactRef, ...] = ()
    checkpoint: ArtifactRef | None = None
    hyperparameters: Mapping[str, Any] = field(default_factory=dict)
    runtime_id: str = "UNRESOLVED"
    runtime_registry_sha256: str = "UNRESOLVED"
    resource_profile: str = "UNRESOLVED"
    code_lock_sha256: str = "UNRESOLVED"
    output_dir: str = ""
    attempt: int = 0
    metadata: Mapping[str, Any] = field(default_factory=dict)

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "campaign_id",
            "task_id",
            "model_id",
            "dataset_ids",
            "split_id",
            "seed",
            "stage",
            "adaptation_regime",
            "fold",
            "action",
            "adapter_command",
            "immutable_inputs",
            "inputs",
            "checkpoint",
            "hyperparameters",
            "runtime_id",
            "runtime_registry_sha256",
            "resource_profile",
            "code_lock_sha256",
            "output_dir",
            "attempt",
            "metadata",
        }
    )
    _REQUIRED: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "campaign_id",
            "task_id",
            "model_id",
            "dataset_ids",
            "split_id",
            "seed",
            "stage",
            "adaptation_regime",
            "fold",
            "action",
            "adapter_command",
            "immutable_inputs",
            "inputs",
            "checkpoint",
            "hyperparameters",
            "runtime_id",
            "runtime_registry_sha256",
            "resource_profile",
            "code_lock_sha256",
            "metadata",
        }
    )

    _ADAPTER_ACTIONS: ClassVar[frozenset[str]] = frozenset(
        {"probe", "prepare", "fit", "predict", "export"}
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "campaign_id",
            "task_id",
            "model_id",
            "split_id",
            "stage",
            "adaptation_regime",
            "runtime_id",
            "resource_profile",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"RunSpec.{field_name}"),
            )
        _set_frozen(
            self,
            "dataset_ids",
            _as_string_tuple(
                self.dataset_ids,
                "RunSpec.dataset_ids",
                allow_empty=self.stage == "terminal_reporting",
            ),
        )
        if not self.dataset_ids and self.stage != "terminal_reporting":
            raise ContractError(
                "only terminal_reporting RunSpecs may omit raw dataset identifiers"
            )
        _set_frozen(self, "seed", _as_int(self.seed, "RunSpec.seed", minimum=0))
        _set_frozen(self, "fold", _as_int(self.fold, "RunSpec.fold", minimum=0))
        actions = _as_string_tuple(
            self.action, "RunSpec.action", allow_empty=False
        )
        unknown_actions = sorted(set(actions).difference(self._ADAPTER_ACTIONS))
        if unknown_actions:
            raise ContractError(
                "RunSpec.action has unsupported adapter action(s): "
                + ", ".join(unknown_actions)
            )
        _set_frozen(self, "action", actions)
        _set_frozen(
            self,
            "adapter_command",
            _as_string_tuple(
                self.adapter_command,
                "RunSpec.adapter_command",
                allow_empty=False,
            ),
        )
        if not isinstance(self.immutable_inputs, Mapping):
            raise ContractError("RunSpec.immutable_inputs must be an object")
        immutable_inputs = {
            _as_string(dataset_id, "RunSpec.immutable_inputs key"): _as_sha256(
                digest,
                f"RunSpec.immutable_inputs[{dataset_id!r}]",
            )
            for dataset_id, digest in self.immutable_inputs.items()
        }
        if set(immutable_inputs) != set(self.dataset_ids):
            raise ContractError(
                "RunSpec.immutable_inputs must cover every and only dataset_id"
            )
        _set_frozen(
            self,
            "immutable_inputs",
            _freeze_json(immutable_inputs, "RunSpec.immutable_inputs"),
        )
        if not isinstance(self.inputs, (tuple, list)):
            raise ContractError("RunSpec.inputs must be an array")
        inputs = tuple(self.inputs)
        if any(not isinstance(item, ArtifactRef) for item in inputs):
            raise ContractError("RunSpec.inputs entries must be ArtifactRef")
        input_paths = [item.path for item in inputs]
        input_roles = [item.role for item in inputs]
        if len(set(input_paths)) != len(input_paths):
            raise ContractError("RunSpec.inputs must not repeat an artifact path")
        if any(role is None for role in input_roles) or len(set(input_roles)) != len(
            input_roles
        ):
            raise ContractError("RunSpec.inputs require unique non-null roles")
        _set_frozen(self, "inputs", inputs)
        if self.checkpoint is not None and not isinstance(self.checkpoint, ArtifactRef):
            raise ContractError("RunSpec.checkpoint must be ArtifactRef or null")
        _set_frozen(
            self,
            "hyperparameters",
            _as_metadata(self.hyperparameters, "RunSpec.hyperparameters"),
        )
        _set_frozen(
            self,
            "runtime_registry_sha256",
            _as_sha256(
                self.runtime_registry_sha256,
                "RunSpec.runtime_registry_sha256",
            ),
        )
        _set_frozen(
            self,
            "code_lock_sha256",
            _as_sha256(self.code_lock_sha256, "RunSpec.code_lock_sha256"),
        )
        _set_frozen(
            self,
            "output_dir",
            _as_string(self.output_dir, "RunSpec.output_dir", allow_empty=True),
        )
        _set_frozen(self, "attempt", _as_int(self.attempt, "RunSpec.attempt", minimum=0))
        _set_frozen(self, "metadata", _as_metadata(self.metadata, "RunSpec.metadata"))

    @property
    def identity_payload(self) -> dict[str, Any]:
        payload = self.to_dict()
        del payload["output_dir"]
        del payload["attempt"]
        return payload

    @property
    def run_id(self) -> str:
        return canonical_sha256(self.identity_payload)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RunSpec":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._REQUIRED,
            label="RunSpec",
        )
        inputs_raw = raw.get("inputs", ())
        if not isinstance(inputs_raw, (list, tuple)):
            raise ContractError("RunSpec.inputs must be an array")
        checkpoint_raw = raw.get("checkpoint")
        return cls(
            schema_version=raw["schema_version"],
            campaign_id=raw["campaign_id"],
            task_id=raw["task_id"],
            model_id=raw["model_id"],
            dataset_ids=raw["dataset_ids"],
            split_id=raw["split_id"],
            seed=raw["seed"],
            stage=raw.get("stage", "benchmark"),
            adaptation_regime=raw.get("adaptation_regime", "common_lane"),
            fold=raw.get("fold", 0),
            action=raw.get("action", ()),
            adapter_command=raw.get("adapter_command", ()),
            immutable_inputs=raw.get("immutable_inputs", {}),
            inputs=tuple(ArtifactRef.from_dict(item) for item in inputs_raw),
            checkpoint=(
                ArtifactRef.from_dict(checkpoint_raw)
                if checkpoint_raw is not None
                else None
            ),
            hyperparameters=raw.get("hyperparameters", {}),
            runtime_id=raw.get("runtime_id", "UNRESOLVED"),
            runtime_registry_sha256=raw.get(
                "runtime_registry_sha256", "UNRESOLVED"
            ),
            resource_profile=raw.get("resource_profile", "UNRESOLVED"),
            code_lock_sha256=raw.get("code_lock_sha256", "UNRESOLVED"),
            output_dir=raw.get("output_dir", ""),
            attempt=raw.get("attempt", 0),
            metadata=raw.get("metadata", {}),
        )


@dataclass(frozen=True, slots=True)
class RunExecutionReceipt(StrictContract):
    """Execution-only evidence that every reviewed adapter action completed."""

    schema_version: str
    candidate_path: str
    campaign_id: str
    plan_sha256: str
    candidate_manifest_sha256: str
    run_id: str
    attempt: int
    wave: str
    task_id: str
    model_id: str
    seed: int
    fold: int
    adaptation_regime: str
    runtime_id: str
    runtime_registry_sha256: str
    adapter_actions: tuple[str, ...]
    action_receipts: tuple[Mapping[str, Any], ...]
    resource_firewall_before_sha256: str | None
    resource_firewall_after_sha256: str | None
    status: RunState
    failure: str | None
    retry_class: str | None

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "candidate_path",
            "campaign_id",
            "plan_sha256",
            "candidate_manifest_sha256",
            "run_id",
            "attempt",
            "wave",
            "task_id",
            "model_id",
            "seed",
            "fold",
            "adaptation_regime",
            "runtime_id",
            "runtime_registry_sha256",
            "adapter_actions",
            "action_receipts",
            "resource_firewall_before_sha256",
            "resource_firewall_after_sha256",
            "status",
            "failure",
            "retry_class",
        }
    )
    _ACTION_RECEIPT_FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "action",
            "status",
            "output_path",
            "adapter_receipt_sha256",
            "output_manifest_sha256",
        }
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "candidate_path",
            "campaign_id",
            "wave",
            "task_id",
            "model_id",
            "adaptation_regime",
            "runtime_id",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(
                    getattr(self, field_name),
                    f"RunExecutionReceipt.{field_name}",
                ),
            )
        if self.schema_version != "masld-bench-run-execution-receipt-v1":
            raise ContractError(
                "RunExecutionReceipt.schema_version must be "
                "masld-bench-run-execution-receipt-v1"
            )
        if not Path(self.candidate_path).is_absolute():
            raise ContractError(
                "RunExecutionReceipt.candidate_path must be an absolute path"
            )
        for field_name in (
            "plan_sha256",
            "candidate_manifest_sha256",
            "run_id",
            "runtime_registry_sha256",
        ):
            _set_frozen(
                self,
                field_name,
                _as_sha256(
                    getattr(self, field_name),
                    f"RunExecutionReceipt.{field_name}",
                ),
            )
        _set_frozen(
            self,
            "attempt",
            _as_int(self.attempt, "RunExecutionReceipt.attempt", minimum=1),
        )
        _set_frozen(
            self,
            "seed",
            _as_int(self.seed, "RunExecutionReceipt.seed", minimum=0),
        )
        _set_frozen(
            self,
            "fold",
            _as_int(self.fold, "RunExecutionReceipt.fold", minimum=0),
        )
        actions = _as_string_tuple(
            self.adapter_actions,
            "RunExecutionReceipt.adapter_actions",
            allow_empty=False,
        )
        if any(action not in RunSpec._ADAPTER_ACTIONS for action in actions):
            raise ContractError("RunExecutionReceipt contains an unsupported adapter action")
        _set_frozen(self, "adapter_actions", actions)
        if not isinstance(self.action_receipts, (list, tuple)):
            raise ContractError("RunExecutionReceipt.action_receipts must be an array")
        action_receipts: list[Mapping[str, Any]] = []
        for index, value in enumerate(self.action_receipts):
            raw = _strict_mapping(
                value,
                allowed=self._ACTION_RECEIPT_FIELDS,
                required=self._ACTION_RECEIPT_FIELDS,
                label=f"RunExecutionReceipt.action_receipts[{index}]",
            )
            action = _as_string(
                raw["action"],
                f"RunExecutionReceipt.action_receipts[{index}].action",
            )
            status = _as_string(
                raw["status"],
                f"RunExecutionReceipt.action_receipts[{index}].status",
            )
            if status != "complete":
                raise ContractError(
                    "RunExecutionReceipt may bind only complete adapter action receipts"
                )
            action_receipts.append(
                MappingProxyType(
                    {
                        "action": action,
                        "status": status,
                        "output_path": _as_string(
                            raw["output_path"],
                            "RunExecutionReceipt action output path",
                        ),
                        "adapter_receipt_sha256": _as_sha256(
                            raw["adapter_receipt_sha256"],
                            "RunExecutionReceipt adapter receipt SHA-256",
                        ),
                        "output_manifest_sha256": _as_sha256(
                            raw["output_manifest_sha256"],
                            "RunExecutionReceipt output manifest SHA-256",
                        ),
                    }
                )
            )
        completed_actions = tuple(item["action"] for item in action_receipts)
        if completed_actions != actions[: len(completed_actions)]:
            raise ContractError(
                "RunExecutionReceipt action receipts must be an ordered action prefix"
            )
        output_paths = [str(item["output_path"]) for item in action_receipts]
        if any(
            Path(path).is_absolute() or ".." in Path(path).parts
            for path in output_paths
        ) or len(set(output_paths)) != len(output_paths):
            raise ContractError(
                "RunExecutionReceipt action output paths must be unique safe relative paths"
            )
        _set_frozen(self, "action_receipts", tuple(action_receipts))
        for field_name in (
            "resource_firewall_before_sha256",
            "resource_firewall_after_sha256",
        ):
            value = getattr(self, field_name)
            _set_frozen(
                self,
                field_name,
                None
                if value is None
                else _as_sha256(value, f"RunExecutionReceipt.{field_name}"),
            )
        _set_frozen(
            self,
            "status",
            _as_enum(self.status, RunState, "RunExecutionReceipt.status"),
        )
        _set_frozen(
            self,
            "failure",
            _optional_string(self.failure, "RunExecutionReceipt.failure"),
        )
        _set_frozen(
            self,
            "retry_class",
            _optional_string(self.retry_class, "RunExecutionReceipt.retry_class"),
        )
        if self.status is RunState.SUCCEEDED:
            if completed_actions != actions:
                raise ContractError(
                    "succeeded RunExecutionReceipt requires every adapter action receipt"
                )
            if self.failure is not None or self.retry_class is not None:
                raise ContractError(
                    "succeeded RunExecutionReceipt may not carry failure metadata"
                )
            if (
                self.resource_firewall_before_sha256 is None
                or self.resource_firewall_after_sha256
                != self.resource_firewall_before_sha256
            ):
                raise ContractError(
                    "succeeded RunExecutionReceipt requires matching pre/post Resource firewalls"
                )
        elif self.failure is None:
            raise ContractError("failed RunExecutionReceipt requires a failure message")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RunExecutionReceipt":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS,
            label="RunExecutionReceipt",
        )
        return cls(**raw)


@dataclass(frozen=True, slots=True)
class PredictionBundle(StrictContract):
    schema_version: str
    bundle_id: str
    run_id: str
    task_id: str
    model_id: str
    dataset_ids: tuple[str, ...]
    split_id: str
    artifacts: tuple[ArtifactRef, ...]
    standardized_table: ArtifactRef
    row_ids: ArtifactRef
    n_predictions: int
    row_id_field: str
    unit_id_field: str
    unit_id_namespace: str
    biological_unit: str
    table_schema_sha256: str
    source_join_key_sha256: str
    format_version: str
    missing_state: MissingState
    metadata: Mapping[str, Any] = field(default_factory=dict)

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "bundle_id",
            "run_id",
            "task_id",
            "model_id",
            "dataset_ids",
            "split_id",
            "artifacts",
            "standardized_table",
            "row_ids",
            "n_predictions",
            "row_id_field",
            "unit_id_field",
            "unit_id_namespace",
            "biological_unit",
            "table_schema_sha256",
            "source_join_key_sha256",
            "format_version",
            "missing_state",
            "metadata",
        }
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "bundle_id",
            "task_id",
            "model_id",
            "split_id",
            "row_id_field",
            "unit_id_field",
            "unit_id_namespace",
            "biological_unit",
            "format_version",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"PredictionBundle.{field_name}"),
            )
        if self.schema_version != "masld-bench-prediction-bundle-v1":
            raise ContractError(
                "PredictionBundle.schema_version must be "
                "masld-bench-prediction-bundle-v1"
            )
        _set_frozen(self, "run_id", _as_sha256(self.run_id, "PredictionBundle.run_id"))
        _set_frozen(
            self,
            "dataset_ids",
            _as_string_tuple(
                self.dataset_ids, "PredictionBundle.dataset_ids", allow_empty=False
            ),
        )
        if not isinstance(self.artifacts, (list, tuple)):
            raise ContractError("PredictionBundle.artifacts must be an array")
        artifacts = tuple(self.artifacts)
        if not artifacts or any(not isinstance(item, ArtifactRef) for item in artifacts):
            raise ContractError("PredictionBundle.artifacts must contain ArtifactRef entries")
        if len({item.path for item in artifacts}) != len(artifacts):
            raise ContractError("PredictionBundle.artifacts must have unique paths")
        if any(item.role is None for item in artifacts) or len(
            {item.role for item in artifacts}
        ) != len(artifacts):
            raise ContractError("PredictionBundle.artifacts require unique non-null roles")
        _set_frozen(self, "artifacts", artifacts)
        for field_name, expected_role in (
            ("standardized_table", f"standardized_prediction_table:{self.task_id}"),
            ("row_ids", f"prediction_row_ids:{self.task_id}"),
        ):
            artifact = getattr(self, field_name)
            if not isinstance(artifact, ArtifactRef):
                raise ContractError(f"PredictionBundle.{field_name} must be ArtifactRef")
            if artifact.role != expected_role:
                raise ContractError(
                    f"PredictionBundle.{field_name} role must be {expected_role!r}"
                )
            matches = [item for item in artifacts if item == artifact]
            if len(matches) != 1:
                raise ContractError(
                    f"PredictionBundle.{field_name} must appear exactly once in artifacts"
                )
        if self.standardized_table.path == self.row_ids.path:
            raise ContractError(
                "PredictionBundle standardized table and row-ID inventory must be distinct"
            )
        _set_frozen(
            self,
            "n_predictions",
            _as_int(self.n_predictions, "PredictionBundle.n_predictions", minimum=1),
        )
        _set_frozen(
            self,
            "table_schema_sha256",
            _as_sha256(
                self.table_schema_sha256,
                "PredictionBundle.table_schema_sha256",
            ),
        )
        expected_join = canonical_sha256(
            {
                "task_id": self.task_id,
                "dataset_ids": list(self.dataset_ids),
                "split_id": self.split_id,
                "row_id_field": self.row_id_field,
                "unit_id_field": self.unit_id_field,
                "unit_id_namespace": self.unit_id_namespace,
                "biological_unit": self.biological_unit,
            }
        )
        _set_frozen(
            self,
            "source_join_key_sha256",
            _as_sha256(
                self.source_join_key_sha256,
                "PredictionBundle.source_join_key_sha256",
            ),
        )
        if self.source_join_key_sha256 != expected_join:
            raise ContractError(
                "PredictionBundle.source_join_key_sha256 does not match its "
                "task/dataset/split/row/unit identity"
            )
        _set_frozen(
            self,
            "missing_state",
            _as_enum(
                self.missing_state, MissingState, "PredictionBundle.missing_state"
            ),
        )
        _set_frozen(
            self, "metadata", _as_metadata(self.metadata, "PredictionBundle.metadata")
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PredictionBundle":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS - {"metadata"},
            label="PredictionBundle",
        )
        artifacts_raw = raw["artifacts"]
        if not isinstance(artifacts_raw, (list, tuple)):
            raise ContractError("PredictionBundle.artifacts must be an array")
        return cls(
            schema_version=raw["schema_version"],
            bundle_id=raw["bundle_id"],
            run_id=raw["run_id"],
            task_id=raw["task_id"],
            model_id=raw["model_id"],
            dataset_ids=raw["dataset_ids"],
            split_id=raw["split_id"],
            artifacts=tuple(ArtifactRef.from_dict(item) for item in artifacts_raw),
            standardized_table=ArtifactRef.from_dict(raw["standardized_table"]),
            row_ids=ArtifactRef.from_dict(raw["row_ids"]),
            n_predictions=raw["n_predictions"],
            row_id_field=raw["row_id_field"],
            unit_id_field=raw["unit_id_field"],
            unit_id_namespace=raw["unit_id_namespace"],
            biological_unit=raw["biological_unit"],
            table_schema_sha256=raw["table_schema_sha256"],
            source_join_key_sha256=raw["source_join_key_sha256"],
            format_version=raw["format_version"],
            missing_state=raw["missing_state"],
            metadata=raw.get("metadata", {}),
        )

    def validate_artifacts(self, root: str | Path) -> tuple[Path, ...]:
        return tuple(
            artifact.validate(root, require_relative=True) for artifact in self.artifacts
        )


@dataclass(frozen=True, slots=True)
class SelectionLock(StrictContract):
    schema_version: str
    lock_id: str
    campaign_id: str
    plan_sha256: str
    candidate_manifest_sha256: str
    registry_sha256: str
    metrics_sha256: str
    locked: bool
    outcomes_unlocked: bool
    release_state: ReleaseState
    candidate_run_ids: tuple[str, ...]
    selected_run_ids: tuple[str, ...]
    task_decisions: tuple[Mapping[str, Any], ...]
    conditional_model_decision: Mapping[str, Any]
    power_decisions: Mapping[str, Any]
    multiplicity_plan: str
    terminal_policy: str
    metadata: Mapping[str, Any] = field(default_factory=dict)

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "lock_id",
            "campaign_id",
            "plan_sha256",
            "candidate_manifest_sha256",
            "registry_sha256",
            "metrics_sha256",
            "locked",
            "outcomes_unlocked",
            "release_state",
            "candidate_run_ids",
            "selected_run_ids",
            "task_decisions",
            "conditional_model_decision",
            "power_decisions",
            "multiplicity_plan",
            "terminal_policy",
            "metadata",
        }
    )

    def __post_init__(self) -> None:
        for field_name in (
            "schema_version",
            "campaign_id",
            "multiplicity_plan",
            "terminal_policy",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"SelectionLock.{field_name}"),
            )
        _set_frozen(self, "lock_id", _as_sha256(self.lock_id, "SelectionLock.lock_id"))
        for field_name in (
            "plan_sha256",
            "candidate_manifest_sha256",
            "registry_sha256",
            "metrics_sha256",
        ):
            _set_frozen(
                self,
                field_name,
                _as_sha256(getattr(self, field_name), f"SelectionLock.{field_name}"),
            )
        _set_frozen(self, "locked", _as_bool(self.locked, "SelectionLock.locked"))
        _set_frozen(
            self,
            "outcomes_unlocked",
            _as_bool(self.outcomes_unlocked, "SelectionLock.outcomes_unlocked"),
        )
        _set_frozen(
            self,
            "release_state",
            _as_enum(self.release_state, ReleaseState, "SelectionLock.release_state"),
        )
        for field_name in ("candidate_run_ids", "selected_run_ids"):
            run_ids = _as_string_tuple(
                getattr(self, field_name),
                f"SelectionLock.{field_name}",
                allow_empty=(field_name == "selected_run_ids"),
            )
            for index, run_id in enumerate(run_ids):
                _as_sha256(run_id, f"SelectionLock.{field_name}[{index}]")
            _set_frozen(self, field_name, run_ids)
        selected_outside_candidates = set(self.selected_run_ids) - set(
            self.candidate_run_ids
        )
        if selected_outside_candidates:
            raise ContractError("selected_run_ids must be a subset of candidate_run_ids")
        if not isinstance(self.task_decisions, (tuple, list)) or not self.task_decisions:
            raise ContractError("SelectionLock.task_decisions must be a non-empty array")
        task_decisions = tuple(
            _as_metadata(item, f"SelectionLock.task_decisions[{index}]")
            for index, item in enumerate(self.task_decisions)
        )
        task_ids = tuple(str(item.get("task_id", "")) for item in task_decisions)
        if any(not task_id for task_id in task_ids) or len(set(task_ids)) != len(task_ids):
            raise ContractError(
                "SelectionLock.task_decisions must have distinct non-empty task_id values"
            )
        _set_frozen(self, "task_decisions", task_decisions)
        _set_frozen(
            self,
            "conditional_model_decision",
            _as_metadata(
                self.conditional_model_decision,
                "SelectionLock.conditional_model_decision",
            ),
        )
        _set_frozen(
            self,
            "power_decisions",
            _as_metadata(self.power_decisions, "SelectionLock.power_decisions"),
        )
        _set_frozen(
            self, "metadata", _as_metadata(self.metadata, "SelectionLock.metadata")
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SelectionLock":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS - {"metadata"},
            label="SelectionLock",
        )
        return cls(
            schema_version=raw["schema_version"],
            lock_id=raw["lock_id"],
            campaign_id=raw["campaign_id"],
            plan_sha256=raw["plan_sha256"],
            candidate_manifest_sha256=raw["candidate_manifest_sha256"],
            registry_sha256=raw["registry_sha256"],
            metrics_sha256=raw["metrics_sha256"],
            locked=raw["locked"],
            outcomes_unlocked=raw["outcomes_unlocked"],
            release_state=raw["release_state"],
            candidate_run_ids=raw["candidate_run_ids"],
            selected_run_ids=raw["selected_run_ids"],
            task_decisions=raw["task_decisions"],
            conditional_model_decision=raw["conditional_model_decision"],
            power_decisions=raw["power_decisions"],
            multiplicity_plan=raw["multiplicity_plan"],
            terminal_policy=raw["terminal_policy"],
            metadata=raw.get("metadata", {}),
        )

    def require_locked(self, *, for_external_scoring: bool = False) -> "SelectionLock":
        if not self.locked:
            raise ContractError(f"selection lock {self.lock_id} is not locked")
        if self.release_state in {ReleaseState.SUPERSEDED, ReleaseState.REJECTED}:
            raise ContractError(
                f"selection lock {self.lock_id} is {self.release_state.value}"
            )
        if for_external_scoring and self.outcomes_unlocked:
            raise ContractError(
                f"selection lock {self.lock_id} was created after outcomes were unlocked"
            )
        return self


Artifact = ArtifactRef
PairingMode = PairingState
MissingnessState = MissingState
AccessTier = AccessState


__all__ = [
    "AccessState",
    "AccessTier",
    "Artifact",
    "ArtifactRef",
    "ContractError",
    "DatasetActivationContract",
    "DatasetManifest",
    "DatasetViewContract",
    "DatasetProvenance",
    "DatasetRole",
    "DatasetSplit",
    "DatasetStatus",
    "ExposureState",
    "ImplementationType",
    "LocusSplitControls",
    "MissingState",
    "MissingnessState",
    "ModelManifest",
    "ModelExecutionContract",
    "PairingMode",
    "PairingState",
    "PredictionFirstOutcomeContract",
    "PredictionBundle",
    "ReleaseClass",
    "ReleaseState",
    "RunSpec",
    "RunExecutionReceipt",
    "RunState",
    "SelectionLock",
    "SplitSpec",
    "StrictContract",
    "TaskSpec",
    "TaskStatus",
]
