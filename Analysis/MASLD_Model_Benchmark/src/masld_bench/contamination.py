"""Deterministic cohort, donor, output file, and checkpoint contamination audits.

This module is deliberately metadata-only.  It accepts salted SHA-256 donor
fingerprints and supplied output-file signatures, never raw donor vectors,
genotypes, metadata, matrices, or images.  Invalid or incomplete evidence is a
hard error; ambiguous checkpoint exposure is an ineligible decision.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from enum import Enum
import re
import unicodedata
from typing import Iterable, Mapping, Sequence
from uuid import UUID

from .contracts import ExposureState
from .hashing import canonical_sha256, require_sha256


class ContaminationAuditError(ValueError):
    """Raised when contamination evidence is incomplete or contradictory."""


class AliasKind(str, Enum):
    ACCESSION = "accession"
    BIOPROJECT = "bioproject"
    SRA_RUN = "sra_run"
    DOI = "doi"
    DONOR = "donor"
    COHORT_RELEASE = "cohort_release"


class FingerprintKind(str, Enum):
    EXPRESSION = "expression"
    GENOTYPE = "genotype"
    METADATA = "metadata"
    SEX = "sex"


class ArtifactKind(str, Enum):
    MATRIX = "matrix"
    IMAGE = "image"


class ExposureConfidence(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    UNKNOWN = "unknown"


CHAMPION_ELIGIBLE_EXPOSURES = frozenset(
    {ExposureState.CLEAN_DECLARED, ExposureState.TARGET_LABEL_UNEXPOSED}
)
SALTED_FINGERPRINT_SCHEME = "sha256_salted_v1"
_UNKNOWN_CORPUS_VALUES = frozenset(
    {
        "not declared",
        "not reported",
        "undeclared",
        "unknown",
        "unreported",
        "unspecified",
    }
)
_WHITESPACE = re.compile(r"\s+")


def _text(value: object, field_name: str) -> str:
    if not isinstance(value, str):
        raise ContaminationAuditError(f"{field_name} must be a string")
    normalized = unicodedata.normalize("NFKC", value).strip()
    if not normalized or "\x00" in normalized:
        raise ContaminationAuditError(f"{field_name} must be a nonempty NUL-free string")
    return normalized


def _identifier(value: object, field_name: str) -> str:
    return _text(value, field_name)


def _parse_date(value: object, field_name: str) -> str:
    normalized = _text(value, field_name)
    try:
        date.fromisoformat(normalized)
    except ValueError as error:
        raise ContaminationAuditError(
            f"{field_name} must be an ISO 8601 date (YYYY-MM-DD)"
        ) from error
    return normalized


def _sha256(value: object, field_name: str) -> str:
    try:
        return require_sha256(value, field_name=field_name)
    except ValueError as error:
        raise ContaminationAuditError(str(error)) from error


def _strict_mapping(
    value: object,
    *,
    fields: frozenset[str],
    label: str,
) -> dict[str, object]:
    if not isinstance(value, Mapping):
        raise ContaminationAuditError(f"{label} must be an object")
    if any(not isinstance(key, str) for key in value):
        raise ContaminationAuditError(f"{label} keys must be strings")
    unknown = sorted(set(value) - fields)
    missing = sorted(fields - set(value))
    if unknown:
        raise ContaminationAuditError(
            f"{label} has unknown field(s): {', '.join(unknown)}"
        )
    if missing:
        raise ContaminationAuditError(
            f"{label} is missing field(s): {', '.join(missing)}"
        )
    return dict(value)


def _sequence(value: object, field_name: str) -> Sequence[object]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise ContaminationAuditError(f"{field_name} must be an array")
    return value


def _normalized_alias(kind: AliasKind, value: object, field_name: str) -> str:
    normalized = _WHITESPACE.sub(" ", _text(value, field_name))
    if kind in {AliasKind.ACCESSION, AliasKind.BIOPROJECT, AliasKind.SRA_RUN}:
        return normalized.replace(" ", "").upper()
    if kind is AliasKind.DOI:
        lowered = normalized.casefold()
        for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
            if lowered.startswith(prefix):
                lowered = lowered[len(prefix) :].strip()
                break
        if "/" not in lowered:
            raise ContaminationAuditError(f"{field_name} is not a qualified DOI")
        return lowered
    lowered = normalized.casefold()
    if kind in {AliasKind.DONOR, AliasKind.COHORT_RELEASE}:
        namespace, separator, local_value = lowered.partition("::")
        if not separator or not namespace.strip() or not local_value.strip():
            raise ContaminationAuditError(
                f"{field_name} must be namespace-qualified as namespace::value"
            )
        return f"{namespace.strip()}::{local_value.strip()}"
    raise ContaminationAuditError(f"unsupported alias kind: {kind}")


def _alias_tuple(kind: AliasKind, value: object, field_name: str) -> tuple[str, ...]:
    aliases = tuple(
        _normalized_alias(kind, item, f"{field_name}[{index}]")
        for index, item in enumerate(_sequence(value, field_name))
    )
    if len(set(aliases)) != len(aliases):
        raise ContaminationAuditError(
            f"{field_name} contains duplicate aliases after normalization"
        )
    return tuple(sorted(aliases))


def _free_alias(value: object, field_name: str) -> str:
    return _WHITESPACE.sub(" ", _text(value, field_name)).casefold()


def _free_alias_tuple(
    value: object, field_name: str, *, allow_empty: bool = True
) -> tuple[str, ...]:
    aliases = tuple(
        _free_alias(item, f"{field_name}[{index}]")
        for index, item in enumerate(_sequence(value, field_name))
    )
    if not allow_empty and not aliases:
        raise ContaminationAuditError(f"{field_name} must not be empty")
    if len(set(aliases)) != len(aliases):
        raise ContaminationAuditError(
            f"{field_name} contains duplicate aliases after normalization"
        )
    return tuple(sorted(aliases))


def _uuid_tuple(value: object, field_name: str) -> tuple[str, ...]:
    identifiers: list[str] = []
    for index, item in enumerate(_sequence(value, field_name)):
        raw = _text(item, f"{field_name}[{index}]")
        try:
            normalized = str(UUID(raw))
        except ValueError as error:
            raise ContaminationAuditError(
                f"{field_name}[{index}] must be a canonical UUID"
            ) from error
        if raw.casefold() != normalized:
            raise ContaminationAuditError(
                f"{field_name}[{index}] must use canonical lowercase UUID form"
            )
        identifiers.append(normalized)
    if len(set(identifiers)) != len(identifiers):
        raise ContaminationAuditError(f"{field_name} contains duplicate UUIDs")
    return tuple(sorted(identifiers))


@dataclass(frozen=True, slots=True, order=True)
class AliasToken:
    kind: AliasKind
    value: str

    def __post_init__(self) -> None:
        try:
            kind = self.kind if isinstance(self.kind, AliasKind) else AliasKind(self.kind)
        except (TypeError, ValueError) as error:
            raise ContaminationAuditError("AliasToken.kind is invalid") from error
        object.__setattr__(self, "kind", kind)
        object.__setattr__(
            self,
            "value",
            _normalized_alias(kind, self.value, "AliasToken.value"),
        )


@dataclass(frozen=True, slots=True)
class CohortIdentityRecord:
    record_id: str
    accession_aliases: tuple[str, ...]
    bioproject_aliases: tuple[str, ...]
    sra_run_aliases: tuple[str, ...]
    doi_aliases: tuple[str, ...]
    donor_aliases: tuple[str, ...]
    cohort_release_aliases: tuple[str, ...]

    _FIELDS = frozenset(
        {
            "record_id",
            "accession_aliases",
            "bioproject_aliases",
            "sra_run_aliases",
            "doi_aliases",
            "donor_aliases",
            "cohort_release_aliases",
        }
    )

    def __post_init__(self) -> None:
        object.__setattr__(self, "record_id", _identifier(self.record_id, "record_id"))
        for field_name, kind in (
            ("accession_aliases", AliasKind.ACCESSION),
            ("bioproject_aliases", AliasKind.BIOPROJECT),
            ("sra_run_aliases", AliasKind.SRA_RUN),
            ("doi_aliases", AliasKind.DOI),
            ("donor_aliases", AliasKind.DONOR),
            ("cohort_release_aliases", AliasKind.COHORT_RELEASE),
        ):
            object.__setattr__(
                self,
                field_name,
                _alias_tuple(kind, getattr(self, field_name), field_name),
            )
        if not self.alias_tokens:
            raise ContaminationAuditError(
                f"cohort identity record {self.record_id!r} has no identity aliases"
            )

    @property
    def alias_tokens(self) -> tuple[AliasToken, ...]:
        tokens: list[AliasToken] = []
        for kind, aliases in (
            (AliasKind.ACCESSION, self.accession_aliases),
            (AliasKind.BIOPROJECT, self.bioproject_aliases),
            (AliasKind.SRA_RUN, self.sra_run_aliases),
            (AliasKind.DOI, self.doi_aliases),
            (AliasKind.DONOR, self.donor_aliases),
            (AliasKind.COHORT_RELEASE, self.cohort_release_aliases),
        ):
            tokens.extend(AliasToken(kind, value) for value in aliases)
        return tuple(sorted(tokens, key=lambda item: (item.kind.value, item.value)))

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> "CohortIdentityRecord":
        return cls(**_strict_mapping(value, fields=cls._FIELDS, label=cls.__name__))


@dataclass(frozen=True, slots=True)
class CohortFamily:
    cohort_family_id: str
    record_ids: tuple[str, ...]
    aliases: tuple[AliasToken, ...]


class _UnionFind:
    def __init__(self, values: Iterable[str]):
        self._parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self._parent[value]
        if parent != value:
            self._parent[value] = self.find(parent)
        return self._parent[value]

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self._parent[max(left_root, right_root)] = min(left_root, right_root)


def resolve_cohort_families(
    records: Sequence[CohortIdentityRecord],
) -> tuple[CohortFamily, ...]:
    """Resolve transitive cohort aliases into stable, order-independent families."""

    by_id: dict[str, CohortIdentityRecord] = {}
    for record in records:
        if not isinstance(record, CohortIdentityRecord):
            raise ContaminationAuditError(
                "resolve_cohort_families requires CohortIdentityRecord values"
            )
        if record.record_id in by_id:
            raise ContaminationAuditError(
                f"duplicate cohort identity record_id: {record.record_id}"
            )
        by_id[record.record_id] = record
    if not by_id:
        return ()

    union = _UnionFind(by_id)
    owners: dict[AliasToken, str] = {}
    for record_id in sorted(by_id):
        for alias in by_id[record_id].alias_tokens:
            prior = owners.setdefault(alias, record_id)
            union.union(record_id, prior)

    members_by_root: dict[str, list[str]] = defaultdict(list)
    for record_id in sorted(by_id):
        members_by_root[union.find(record_id)].append(record_id)

    families: list[CohortFamily] = []
    seen_family_ids: set[str] = set()
    for member_ids in members_by_root.values():
        members = tuple(sorted(member_ids))
        aliases = tuple(
            sorted(
                {alias for record_id in members for alias in by_id[record_id].alias_tokens},
                key=lambda item: (item.kind.value, item.value),
            )
        )
        family_hash = canonical_sha256(
            {
                "schema": "masld-bench-cohort-family-v1",
                "aliases": [
                    {"kind": alias.kind.value, "value": alias.value} for alias in aliases
                ],
            }
        )
        family_id = f"cohort-family-{family_hash}"
        if family_id in seen_family_ids:
            raise ContaminationAuditError("cohort family hash collision")
        seen_family_ids.add(family_id)
        families.append(CohortFamily(family_id, members, aliases))
    return tuple(sorted(families, key=lambda item: item.cohort_family_id))


@dataclass(frozen=True, slots=True)
class SaltedDonorFingerprintRecord:
    subject_record_id: str
    fingerprint_kind: FingerprintKind
    fingerprint_sha256: str
    salt_commitment_sha256: str
    scheme: str = SALTED_FINGERPRINT_SCHEME

    _FIELDS = frozenset(
        {
            "subject_record_id",
            "fingerprint_kind",
            "fingerprint_sha256",
            "salt_commitment_sha256",
            "scheme",
        }
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "subject_record_id",
            _identifier(self.subject_record_id, "subject_record_id"),
        )
        try:
            kind = (
                self.fingerprint_kind
                if isinstance(self.fingerprint_kind, FingerprintKind)
                else FingerprintKind(self.fingerprint_kind)
            )
        except (TypeError, ValueError) as error:
            raise ContaminationAuditError("fingerprint_kind is invalid") from error
        object.__setattr__(self, "fingerprint_kind", kind)
        object.__setattr__(
            self,
            "fingerprint_sha256",
            _sha256(self.fingerprint_sha256, "fingerprint_sha256"),
        )
        object.__setattr__(
            self,
            "salt_commitment_sha256",
            _sha256(self.salt_commitment_sha256, "salt_commitment_sha256"),
        )
        scheme = _text(self.scheme, "scheme")
        if scheme != SALTED_FINGERPRINT_SCHEME:
            raise ContaminationAuditError(
                f"scheme must be exactly {SALTED_FINGERPRINT_SCHEME!r}"
            )
        object.__setattr__(self, "scheme", scheme)

    @classmethod
    def from_dict(
        cls, value: Mapping[str, object]
    ) -> "SaltedDonorFingerprintRecord":
        return cls(**_strict_mapping(value, fields=cls._FIELDS, label=cls.__name__))


@dataclass(frozen=True, slots=True)
class DonorDuplicateComponent:
    component_id: str
    subject_record_ids: tuple[str, ...]
    matching_fingerprint_kinds: tuple[FingerprintKind, ...]


def donor_duplicate_components(
    records: Sequence[SaltedDonorFingerprintRecord],
) -> tuple[DonorDuplicateComponent, ...]:
    """Group donors sharing salted signatures and reject within-component conflicts."""

    if not records:
        return ()
    for record in records:
        if not isinstance(record, SaltedDonorFingerprintRecord):
            raise ContaminationAuditError(
                "donor_duplicate_components requires SaltedDonorFingerprintRecord values"
            )
    salt_commitments = {record.salt_commitment_sha256 for record in records}
    if len(salt_commitments) != 1:
        raise ContaminationAuditError(
            "all donor fingerprints in one audit must use the same salt commitment"
        )

    by_subject_kind: dict[
        tuple[str, FingerprintKind], SaltedDonorFingerprintRecord
    ] = {}
    for record in records:
        key = (record.subject_record_id, record.fingerprint_kind)
        if key in by_subject_kind:
            raise ContaminationAuditError(
                "duplicate or conflicting fingerprint for "
                f"{record.subject_record_id}/{record.fingerprint_kind.value}"
            )
        by_subject_kind[key] = record

    subjects = sorted({record.subject_record_id for record in records})
    union = _UnionFind(subjects)
    signature_owners: dict[tuple[FingerprintKind, str], str] = {}
    for record in sorted(
        records,
        key=lambda item: (
            item.fingerprint_kind.value,
            item.fingerprint_sha256,
            item.subject_record_id,
        ),
    ):
        if record.fingerprint_kind not in {
            FingerprintKind.EXPRESSION,
            FingerprintKind.GENOTYPE,
        }:
            continue
        signature = (record.fingerprint_kind, record.fingerprint_sha256)
        prior = signature_owners.setdefault(signature, record.subject_record_id)
        if prior != record.subject_record_id:
            union.union(prior, record.subject_record_id)

    members_by_root: dict[str, list[str]] = defaultdict(list)
    for subject in subjects:
        members_by_root[union.find(subject)].append(subject)

    components: list[DonorDuplicateComponent] = []
    for members_raw in members_by_root.values():
        members = tuple(sorted(members_raw))
        if len(members) < 2:
            continue
        member_set = set(members)
        signature_subjects_by_kind: dict[
            FingerprintKind, dict[str, set[str]]
        ] = defaultdict(lambda: defaultdict(set))
        for (subject, kind), record in by_subject_kind.items():
            if subject in member_set:
                signature_subjects_by_kind[kind][record.fingerprint_sha256].add(subject)
        conflicts = sorted(
            kind.value
            for kind, signature_subjects in signature_subjects_by_kind.items()
            if kind in {FingerprintKind.GENOTYPE, FingerprintKind.SEX}
            and len(signature_subjects) > 1
        )
        if conflicts:
            raise ContaminationAuditError(
                "duplicate donor component has conflicting stable salted fingerprint kind(s): "
                + ", ".join(conflicts)
            )
        matching_kinds = tuple(
            sorted(
                (
                    kind for kind, signature_subjects in signature_subjects_by_kind.items()
                    if any(len(subjects) > 1 for subjects in signature_subjects.values())
                ),
                key=lambda item: item.value,
            )
        )
        component_hash = canonical_sha256(
            {
                "schema": "masld-bench-donor-duplicate-component-v1",
                "subjects": members,
                "salt_commitment_sha256": next(iter(salt_commitments)),
                "salted_fingerprints": [
                    {
                        "subject_record_id": record.subject_record_id,
                        "kind": record.fingerprint_kind.value,
                        "fingerprint_sha256": record.fingerprint_sha256,
                    }
                    for record in sorted(
                        (
                            record
                            for record in records
                            if record.subject_record_id in member_set
                        ),
                        key=lambda item: (
                            item.subject_record_id,
                            item.fingerprint_kind.value,
                        ),
                    )
                ],
                "matching_fingerprint_kinds": [kind.value for kind in matching_kinds],
            }
        )
        components.append(
            DonorDuplicateComponent(
                component_id=f"donor-duplicate-{component_hash}",
                subject_record_ids=members,
                matching_fingerprint_kinds=matching_kinds,
            )
        )
    return tuple(sorted(components, key=lambda item: item.component_id))


@dataclass(frozen=True, slots=True)
class ArtifactSignatureRecord:
    artifact_id: str
    artifact_kind: ArtifactKind
    exact_sha256: str
    supplied_near_signature_sha256: str | None

    _FIELDS = frozenset(
        {
            "artifact_id",
            "artifact_kind",
            "exact_sha256",
            "supplied_near_signature_sha256",
        }
    )

    def __post_init__(self) -> None:
        object.__setattr__(self, "artifact_id", _identifier(self.artifact_id, "artifact_id"))
        try:
            kind = (
                self.artifact_kind
                if isinstance(self.artifact_kind, ArtifactKind)
                else ArtifactKind(self.artifact_kind)
            )
        except (TypeError, ValueError) as error:
            raise ContaminationAuditError("artifact_kind is invalid") from error
        object.__setattr__(self, "artifact_kind", kind)
        object.__setattr__(self, "exact_sha256", _sha256(self.exact_sha256, "exact_sha256"))
        if self.supplied_near_signature_sha256 is not None:
            object.__setattr__(
                self,
                "supplied_near_signature_sha256",
                _sha256(
                    self.supplied_near_signature_sha256,
                    "supplied_near_signature_sha256",
                ),
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> "ArtifactSignatureRecord":
        return cls(**_strict_mapping(value, fields=cls._FIELDS, label=cls.__name__))


@dataclass(frozen=True, slots=True)
class ArtifactDuplicateComponent:
    component_id: str
    artifact_kind: ArtifactKind
    artifact_ids: tuple[str, ...]
    evidence: tuple[str, ...]


def artifact_duplicate_components(
    records: Sequence[ArtifactSignatureRecord],
) -> tuple[ArtifactDuplicateComponent, ...]:
    """Group exact and supplied near-duplicate matrix or image signatures."""

    by_id: dict[str, ArtifactSignatureRecord] = {}
    exact_kinds: dict[str, ArtifactKind] = {}
    for record in records:
        if not isinstance(record, ArtifactSignatureRecord):
            raise ContaminationAuditError(
                "artifact_duplicate_components requires ArtifactSignatureRecord values"
            )
        if record.artifact_id in by_id:
            raise ContaminationAuditError(f"duplicate artifact_id: {record.artifact_id}")
        prior_kind = exact_kinds.setdefault(record.exact_sha256, record.artifact_kind)
        if prior_kind is not record.artifact_kind:
            raise ContaminationAuditError(
                "one exact artifact SHA-256 is classified as both matrix and image"
            )
        by_id[record.artifact_id] = record
    if not by_id:
        return ()

    union = _UnionFind(by_id)
    exact_owner: dict[tuple[ArtifactKind, str], str] = {}
    near_owner: dict[tuple[ArtifactKind, str], str] = {}
    for artifact_id in sorted(by_id):
        record = by_id[artifact_id]
        exact_key = (record.artifact_kind, record.exact_sha256)
        union.union(artifact_id, exact_owner.setdefault(exact_key, artifact_id))
        if record.supplied_near_signature_sha256 is not None:
            near_key = (record.artifact_kind, record.supplied_near_signature_sha256)
            union.union(artifact_id, near_owner.setdefault(near_key, artifact_id))

    members_by_root: dict[str, list[str]] = defaultdict(list)
    for artifact_id in sorted(by_id):
        members_by_root[union.find(artifact_id)].append(artifact_id)

    components: list[ArtifactDuplicateComponent] = []
    for member_ids in members_by_root.values():
        members = tuple(sorted(member_ids))
        if len(members) < 2:
            continue
        kinds = {by_id[member].artifact_kind for member in members}
        if len(kinds) != 1:
            raise ContaminationAuditError("artifact duplicate component crosses artifact kinds")
        kind = next(iter(kinds))
        evidence: set[str] = set()
        exact_counts: dict[str, int] = defaultdict(int)
        near_counts: dict[str, int] = defaultdict(int)
        for member in members:
            record = by_id[member]
            exact_counts[record.exact_sha256] += 1
            if record.supplied_near_signature_sha256 is not None:
                near_counts[record.supplied_near_signature_sha256] += 1
        if any(count > 1 for count in exact_counts.values()):
            evidence.add("exact_sha256")
        if any(count > 1 for count in near_counts.values()):
            evidence.add("supplied_near_signature_sha256")
        component_hash = canonical_sha256(
            {
                "schema": "masld-bench-artifact-duplicate-component-v1",
                "artifact_kind": kind.value,
                "artifacts": members,
                "signatures": [
                    {
                        "artifact_id": member,
                        "exact_sha256": by_id[member].exact_sha256,
                        "supplied_near_signature_sha256": by_id[
                            member
                        ].supplied_near_signature_sha256,
                    }
                    for member in members
                ],
                "evidence": sorted(evidence),
            }
        )
        components.append(
            ArtifactDuplicateComponent(
                component_id=f"artifact-duplicate-{component_hash}",
                artifact_kind=kind,
                artifact_ids=members,
                evidence=tuple(sorted(evidence)),
            )
        )
    return tuple(sorted(components, key=lambda item: item.component_id))


@dataclass(frozen=True, slots=True)
class CheckpointAuditRecord:
    schema_version: str
    model_id: str
    checkpoint_revision: str
    checkpoint_sha256: str
    release_date: str
    training_cutoff: str
    declared_corpora: tuple[str, ...]
    cellxgene_uuids: tuple[str, ...]
    accession_aliases: tuple[str, ...]
    corpus_declaration_complete: bool
    exposure_state: ExposureState
    exposure_confidence: ExposureConfidence

    _FIELDS = frozenset(
        {
            "schema_version",
            "model_id",
            "checkpoint_revision",
            "checkpoint_sha256",
            "release_date",
            "training_cutoff",
            "declared_corpora",
            "cellxgene_uuids",
            "accession_aliases",
            "corpus_declaration_complete",
            "exposure_state",
            "exposure_confidence",
        }
    )

    def __post_init__(self) -> None:
        for field_name in ("schema_version", "model_id", "checkpoint_revision"):
            value = _identifier(getattr(self, field_name), field_name)
            if value in {"UNRESOLVED", "UNKNOWN", "NOT_APPLICABLE"}:
                raise ContaminationAuditError(f"{field_name} must be exact")
            object.__setattr__(self, field_name, value)
        object.__setattr__(
            self,
            "checkpoint_sha256",
            _sha256(self.checkpoint_sha256, "checkpoint_sha256"),
        )
        release_date = _parse_date(self.release_date, "release_date")
        training_cutoff = _parse_date(self.training_cutoff, "training_cutoff")
        if date.fromisoformat(training_cutoff) > date.fromisoformat(release_date):
            raise ContaminationAuditError("training_cutoff must not follow release_date")
        object.__setattr__(self, "release_date", release_date)
        object.__setattr__(self, "training_cutoff", training_cutoff)
        object.__setattr__(
            self,
            "declared_corpora",
            _free_alias_tuple(self.declared_corpora, "declared_corpora", allow_empty=False),
        )
        object.__setattr__(
            self,
            "cellxgene_uuids",
            _uuid_tuple(self.cellxgene_uuids, "cellxgene_uuids"),
        )
        object.__setattr__(
            self,
            "accession_aliases",
            _alias_tuple(AliasKind.ACCESSION, self.accession_aliases, "accession_aliases"),
        )
        if not isinstance(self.corpus_declaration_complete, bool):
            raise ContaminationAuditError("corpus_declaration_complete must be boolean")
        try:
            exposure_state = (
                self.exposure_state
                if isinstance(self.exposure_state, ExposureState)
                else ExposureState(self.exposure_state)
            )
        except (TypeError, ValueError) as error:
            raise ContaminationAuditError("exposure_state is invalid") from error
        object.__setattr__(self, "exposure_state", exposure_state)
        try:
            confidence = (
                self.exposure_confidence
                if isinstance(self.exposure_confidence, ExposureConfidence)
                else ExposureConfidence(self.exposure_confidence)
            )
        except (TypeError, ValueError) as error:
            raise ContaminationAuditError("exposure_confidence is invalid") from error
        object.__setattr__(self, "exposure_confidence", confidence)

    @property
    def audit_id(self) -> str:
        return canonical_sha256(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> "CheckpointAuditRecord":
        return cls(**_strict_mapping(value, fields=cls._FIELDS, label=cls.__name__))


@dataclass(frozen=True, slots=True)
class ScBaseCountSnapshot:
    snapshot_id: str
    release_date: str
    corpus_aliases: tuple[str, ...]
    accession_aliases: tuple[str, ...]
    cellxgene_uuids: tuple[str, ...]

    _FIELDS = frozenset(
        {
            "snapshot_id",
            "release_date",
            "corpus_aliases",
            "accession_aliases",
            "cellxgene_uuids",
        }
    )

    def __post_init__(self) -> None:
        object.__setattr__(self, "snapshot_id", _identifier(self.snapshot_id, "snapshot_id"))
        object.__setattr__(self, "release_date", _parse_date(self.release_date, "release_date"))
        object.__setattr__(
            self,
            "corpus_aliases",
            _free_alias_tuple(self.corpus_aliases, "corpus_aliases", allow_empty=False),
        )
        if not any("scbasecount" in alias for alias in self.corpus_aliases):
            raise ContaminationAuditError(
                "ScBaseCountSnapshot.corpus_aliases must identify scBaseCount"
            )
        object.__setattr__(
            self,
            "accession_aliases",
            _alias_tuple(AliasKind.ACCESSION, self.accession_aliases, "accession_aliases"),
        )
        object.__setattr__(
            self,
            "cellxgene_uuids",
            _uuid_tuple(self.cellxgene_uuids, "cellxgene_uuids"),
        )

    @property
    def snapshot_hash(self) -> str:
        return canonical_sha256(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> "ScBaseCountSnapshot":
        return cls(**_strict_mapping(value, fields=cls._FIELDS, label=cls.__name__))


@dataclass(frozen=True, slots=True)
class CheckpointEligibility:
    decision_id: str
    checkpoint_audit_id: str
    target_accession: str
    eligible: bool
    reasons: tuple[str, ...]
    forbidden_snapshot_ids: tuple[str, ...]


def assess_checkpoint_champion_eligibility(
    record: CheckpointAuditRecord,
    *,
    scbasecount_snapshots: Sequence[ScBaseCountSnapshot],
    target_accession: str = "GSE289173",
) -> CheckpointEligibility:
    """Apply the fail-closed held-back evaluation checkpoint exposure policy."""

    if not isinstance(record, CheckpointAuditRecord):
        raise ContaminationAuditError("record must be a CheckpointAuditRecord")
    target = _normalized_alias(
        AliasKind.ACCESSION, target_accession, "target_accession"
    )
    snapshots_by_id: dict[str, ScBaseCountSnapshot] = {}
    alias_owner: dict[str, str] = {}
    uuid_owner: dict[str, str] = {}
    for snapshot in scbasecount_snapshots:
        if not isinstance(snapshot, ScBaseCountSnapshot):
            raise ContaminationAuditError(
                "scbasecount_snapshots must contain ScBaseCountSnapshot values"
            )
        if snapshot.snapshot_id in snapshots_by_id:
            raise ContaminationAuditError(
                f"duplicate scBaseCount snapshot_id: {snapshot.snapshot_id}"
            )
        snapshots_by_id[snapshot.snapshot_id] = snapshot
        for alias in snapshot.corpus_aliases:
            prior = alias_owner.setdefault(alias, snapshot.snapshot_id)
            if prior != snapshot.snapshot_id:
                raise ContaminationAuditError(
                    f"scBaseCount corpus alias maps to multiple snapshots: {alias}"
                )
        for identifier in snapshot.cellxgene_uuids:
            prior = uuid_owner.setdefault(identifier, snapshot.snapshot_id)
            if prior != snapshot.snapshot_id:
                raise ContaminationAuditError(
                    "CELLxGENE UUID maps to multiple scBaseCount snapshots: "
                    + identifier
                )

    target_snapshots = tuple(
        sorted(
            (
                snapshot
                for snapshot in snapshots_by_id.values()
                if date.fromisoformat(snapshot.release_date).year == 2026
                and target in snapshot.accession_aliases
            ),
            key=lambda item: (item.release_date, item.snapshot_id),
        )
    )
    reasons: set[str] = set()
    forbidden_ids: set[str] = set()
    if not target_snapshots:
        reasons.add("target_2026_scbasecount_snapshot_registry_missing")
        anchor_date: date | None = None
    else:
        anchor_date = date.fromisoformat(target_snapshots[0].release_date)

    if record.exposure_state not in CHAMPION_ELIGIBLE_EXPOSURES:
        reasons.add(
            "exposure_state_not_champion_eligible:" + record.exposure_state.value
        )
    if record.exposure_confidence is not ExposureConfidence.HIGH:
        reasons.add(
            "exposure_confidence_not_high:" + record.exposure_confidence.value
        )
    if not record.corpus_declaration_complete:
        reasons.add("corpus_declaration_incomplete")
    for corpus in record.declared_corpora:
        if corpus in _UNKNOWN_CORPUS_VALUES:
            reasons.add("unknown_or_undeclared_corpus:" + corpus.replace(" ", "_"))

    if target in record.accession_aliases:
        reasons.add("target_accession_declared:" + target)

    matched_snapshot_ids: set[str] = set()
    for corpus in record.declared_corpora:
        snapshot_id = alias_owner.get(corpus)
        if snapshot_id is not None:
            matched_snapshot_ids.add(snapshot_id)
        elif "scbasecount" in corpus:
            reasons.add("unregistered_scbasecount_snapshot:" + corpus)
    for identifier in record.cellxgene_uuids:
        snapshot_id = uuid_owner.get(identifier)
        if snapshot_id is not None:
            matched_snapshot_ids.add(snapshot_id)

    declares_cellxgene = any("cellxgene" in corpus for corpus in record.declared_corpora)
    if declares_cellxgene and not record.cellxgene_uuids:
        reasons.add("cellxgene_corpus_without_exact_uuids")

    if anchor_date is not None:
        for snapshot_id in matched_snapshot_ids:
            snapshot = snapshots_by_id[snapshot_id]
            if date.fromisoformat(snapshot.release_date) >= anchor_date:
                forbidden_ids.add(snapshot_id)
                reasons.add("forbidden_scbasecount_snapshot:" + snapshot_id)
            if date.fromisoformat(record.training_cutoff) < date.fromisoformat(
                snapshot.release_date
            ):
                reasons.add("snapshot_after_declared_training_cutoff:" + snapshot_id)

    ordered_reasons = tuple(sorted(reasons))
    ordered_forbidden = tuple(sorted(forbidden_ids))
    decision_payload = {
        "schema": "masld-bench-checkpoint-eligibility-v1",
        "checkpoint_audit_id": record.audit_id,
        "target_accession": target,
        "snapshot_hashes": [
            snapshots_by_id[snapshot_id].snapshot_hash
            for snapshot_id in sorted(snapshots_by_id)
        ],
        "eligible": not ordered_reasons,
        "reasons": ordered_reasons,
        "forbidden_snapshot_ids": ordered_forbidden,
    }
    return CheckpointEligibility(
        decision_id=canonical_sha256(decision_payload),
        checkpoint_audit_id=record.audit_id,
        target_accession=target,
        eligible=not ordered_reasons,
        reasons=ordered_reasons,
        forbidden_snapshot_ids=ordered_forbidden,
    )


__all__ = [
    "AliasKind",
    "AliasToken",
    "ArtifactDuplicateComponent",
    "ArtifactKind",
    "ArtifactSignatureRecord",
    "CHAMPION_ELIGIBLE_EXPOSURES",
    "CheckpointAuditRecord",
    "CheckpointEligibility",
    "CohortFamily",
    "CohortIdentityRecord",
    "ContaminationAuditError",
    "DonorDuplicateComponent",
    "ExposureConfidence",
    "FingerprintKind",
    "SALTED_FINGERPRINT_SCHEME",
    "SaltedDonorFingerprintRecord",
    "ScBaseCountSnapshot",
    "artifact_duplicate_components",
    "assess_checkpoint_champion_eligibility",
    "donor_duplicate_components",
    "resolve_cohort_families",
]
