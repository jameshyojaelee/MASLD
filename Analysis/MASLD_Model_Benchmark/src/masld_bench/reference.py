"""Fail-closed genome and annotation reference contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
import gzip
from pathlib import Path
from typing import Any, ClassVar, Mapping

from .contracts import (
    ContractError,
    StrictContract,
    _as_sha256,
    _as_string,
    _as_string_tuple,
    _set_frozen,
    _strict_mapping,
)
from .hashing import sha256_file


class ReferenceError(ContractError):
    """Raised when a reference bundle is unsafe, corrupt, or mismatched."""


KNOWN_GENCODE_V49_ANNOTATION_SHA256 = (
    "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
)
KNOWN_GRCH38_P14_GENOME_SHA256 = (
    "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215"
)
KNOWN_TRUNCATED_GENCODE_V49_PATHS = frozenset(
    {
        "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/GRCh38.p14.genome.fa.gz",
        "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.transcripts.fa.gz",
    }
)
KNOWN_TRUNCATED_GENCODE_V49_SHA256 = frozenset(
    {
        "925a7184ff41a6dcaac15db88e953256a1c40576fc4fb42dae1958630bd873e2",
        "a4456e191e0b3d1f6e55c33f8482639b8b2bd02dfbc99af2f4fa19e96b65a6be",
    }
)


def _resolve_reference_path(path: str, root: str | Path | None) -> Path:
    configured = Path(path)
    if configured.is_absolute():
        candidate = configured
    elif root is not None:
        candidate = Path(root) / configured
    else:
        candidate = configured
    try:
        resolved = candidate.resolve(strict=True)
    except OSError as error:
        raise ReferenceError(f"reference is missing or unreadable: {candidate}") from error
    if not resolved.is_file():
        raise ReferenceError(f"reference is not a regular file: {resolved}")
    return resolved


def _reject_known_truncated(path: Path, declared_sha256: str) -> None:
    if (
        path.as_posix() in KNOWN_TRUNCATED_GENCODE_V49_PATHS
        or declared_sha256 in KNOWN_TRUNCATED_GENCODE_V49_SHA256
    ):
        raise ReferenceError(f"reference is a known truncated GENCODE v49 file: {path}")


def _validate_integrity_label(value: str, *, gzip_required: bool) -> None:
    normalized = value.casefold().replace("-", "_")
    if normalized in {"unknown", "unchecked", "unresolved", "none"}:
        raise ReferenceError(f"reference integrity is not verified: {value!r}")
    if "sha" not in normalized:
        raise ReferenceError("reference integrity declaration must include SHA-256")
    if gzip_required and ("gzip" not in normalized and "gz" not in normalized):
        raise ReferenceError("gzip reference integrity declaration must include gzip")


def _stream_validate_gzip(path: Path, *, chunk_size: int = 8 * 1024 * 1024) -> None:
    decompressed_bytes = 0
    try:
        with gzip.open(path, "rb") as handle:
            while chunk := handle.read(chunk_size):
                decompressed_bytes += len(chunk)
    except (gzip.BadGzipFile, EOFError, OSError) as error:
        raise ReferenceError(f"gzip integrity check failed for {path}: {error}") from error
    if decompressed_bytes == 0:
        raise ReferenceError(f"gzip reference decompresses to an empty file: {path}")


def _fasta_contigs(path: Path) -> frozenset[str]:
    contigs: set[str] = set()
    try:
        with gzip.open(path, "rt", encoding="ascii") as handle:
            for line in handle:
                if line.startswith(">"):
                    identifier = line[1:].split(maxsplit=1)[0].strip()
                    if not identifier or identifier in contigs:
                        raise ReferenceError(
                            f"FASTA has an empty or duplicate contig header: {path}"
                        )
                    contigs.add(identifier)
    except (gzip.BadGzipFile, EOFError, OSError, UnicodeError) as error:
        raise ReferenceError(f"cannot read FASTA contigs from {path}: {error}") from error
    if not contigs:
        raise ReferenceError(f"FASTA has no contig headers: {path}")
    return frozenset(contigs)


def _annotation_contigs(path: Path) -> frozenset[str]:
    contigs: set[str] = set()
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t", maxsplit=1)
                if len(fields) != 2 or not fields[0]:
                    raise ReferenceError(f"annotation has a malformed feature row: {path}")
                contigs.add(fields[0])
    except (gzip.BadGzipFile, EOFError, OSError, UnicodeError) as error:
        raise ReferenceError(f"cannot read annotation contigs from {path}: {error}") from error
    if not contigs:
        raise ReferenceError(f"annotation has no feature contigs: {path}")
    return frozenset(contigs)


def _validate_reference_file(
    *,
    path: str,
    declared_sha256: str,
    integrity: str,
    gzip_required: bool,
    root: str | Path | None,
) -> Path:
    source = _resolve_reference_path(path, root)
    _reject_known_truncated(source, declared_sha256)
    _validate_integrity_label(integrity, gzip_required=gzip_required)
    actual_sha256 = sha256_file(source)
    if actual_sha256 != declared_sha256:
        raise ReferenceError(
            f"reference SHA-256 mismatch for {source}: expected {declared_sha256}, got {actual_sha256}"
        )
    if gzip_required:
        _stream_validate_gzip(source)
    return source


@dataclass(frozen=True, slots=True)
class GenomeReference(StrictContract):
    path: str
    sha256: str
    compression: str
    integrity: str

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {"path", "sha256", "compression", "integrity"}
    )

    def __post_init__(self) -> None:
        for field_name in ("path", "compression", "integrity"):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"GenomeReference.{field_name}"),
            )
        _set_frozen(self, "sha256", _as_sha256(self.sha256, "GenomeReference.sha256"))
        if self.compression not in {"gzip", "none"}:
            raise ReferenceError("GenomeReference.compression must be 'gzip' or 'none'")
        if self.compression == "gzip" and not self.path.endswith(".gz"):
            raise ReferenceError("gzip genome reference path must end in .gz")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "GenomeReference":
        raw = _strict_mapping(
            value, allowed=cls._FIELDS, required=cls._FIELDS, label="GenomeReference"
        )
        return cls(**raw)

    def validate(self, root: str | Path | None = None) -> Path:
        return _validate_reference_file(
            path=self.path,
            declared_sha256=self.sha256,
            integrity=self.integrity,
            gzip_required=self.compression == "gzip",
            root=root,
        )


@dataclass(frozen=True, slots=True)
class AnnotationReference(StrictContract):
    path: str
    sha256: str
    format: str
    integrity: str

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {"path", "sha256", "format", "integrity"}
    )

    def __post_init__(self) -> None:
        for field_name in ("path", "format", "integrity"):
            _set_frozen(
                self,
                field_name,
                _as_string(
                    getattr(self, field_name), f"AnnotationReference.{field_name}"
                ),
            )
        _set_frozen(
            self, "sha256", _as_sha256(self.sha256, "AnnotationReference.sha256")
        )
        if self.format not in {"gtf", "gff3", "gtf_gzip", "gff3_gzip"}:
            raise ReferenceError(
                "AnnotationReference.format must describe gzip GTF or GFF3"
            )
        if not self.path.endswith(".gz"):
            raise ReferenceError("annotation reference must be gzip-compressed (.gz)")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "AnnotationReference":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS,
            label="AnnotationReference",
        )
        return cls(**raw)

    def validate(self, root: str | Path | None = None) -> Path:
        return _validate_reference_file(
            path=self.path,
            declared_sha256=self.sha256,
            integrity=self.integrity,
            gzip_required=True,
            root=root,
        )


@dataclass(frozen=True, slots=True)
class CoordinateContract(StrictContract):
    bed_system: str
    contig_policy: str
    analysis_contigs: tuple[str, ...]
    excluded_contig_policy: str
    variant_normalization: str
    liftover_policy: str
    native_ld_policy: str

    _FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "bed_system",
            "contig_policy",
            "analysis_contigs",
            "excluded_contig_policy",
            "variant_normalization",
            "liftover_policy",
            "native_ld_policy",
        }
    )

    def __post_init__(self) -> None:
        for field_name in self._FIELDS - {"analysis_contigs"}:
            _set_frozen(
                self,
                field_name,
                _as_string(
                    getattr(self, field_name), f"CoordinateContract.{field_name}"
                ),
            )
        _set_frozen(
            self,
            "analysis_contigs",
            _as_string_tuple(
                self.analysis_contigs,
                "CoordinateContract.analysis_contigs",
                allow_empty=False,
            ),
        )
        if self.bed_system != "zero_based_half_open":
            raise ReferenceError(
                "CoordinateContract.bed_system must be zero_based_half_open"
            )
        if self.excluded_contig_policy != "reject":
            raise ReferenceError(
                "CoordinateContract.excluded_contig_policy must be reject"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "CoordinateContract":
        raw = _strict_mapping(
            value,
            allowed=cls._FIELDS,
            required=cls._FIELDS,
            label="CoordinateContract",
        )
        return cls(**raw)

    def require_analysis_contig(self, contig: str) -> str:
        normalized = _as_string(contig, "contig")
        if normalized not in self.analysis_contigs:
            raise ReferenceError(
                f"contig is excluded from the analysis coordinate contract: {normalized}"
            )
        return normalized


@dataclass(frozen=True, slots=True)
class ReferenceBundle(StrictContract):
    bundle_id: str
    assembly: str
    assembly_patch: str
    annotation_release: str
    genome: GenomeReference
    annotation: AnnotationReference
    sequence_extraction_ready: bool
    sequence_extraction_blockers: tuple[str, ...]
    coordinate_contract: CoordinateContract
    prohibited_paths: tuple[str, ...]
    prohibited_reason: str
    _base_dir: Path | None = field(default=None, repr=False, compare=False)

    _ROOT_FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "schema_version",
            "reference",
            "coordinate_contract",
            "semantics",
            "prohibited_references",
            "redistribution",
            "cluster",
            "profiles",
        }
    )
    _ROOT_REQUIRED: ClassVar[frozenset[str]] = frozenset(
        {"reference", "coordinate_contract", "prohibited_references"}
    )
    _REFERENCE_FIELDS: ClassVar[frozenset[str]] = frozenset(
        {
            "bundle_id",
            "assembly",
            "assembly_patch",
            "annotation_release",
            "genome",
            "annotation",
            "sequence_extraction_ready",
            "sequence_extraction_blockers",
        }
    )
    _PROHIBITED_FIELDS: ClassVar[frozenset[str]] = frozenset({"paths", "reason"})

    def __post_init__(self) -> None:
        for field_name in (
            "bundle_id",
            "assembly",
            "assembly_patch",
            "annotation_release",
        ):
            _set_frozen(
                self,
                field_name,
                _as_string(getattr(self, field_name), f"ReferenceBundle.{field_name}"),
            )
        if not isinstance(self.genome, GenomeReference):
            raise ReferenceError("ReferenceBundle.genome must be GenomeReference")
        if not isinstance(self.annotation, AnnotationReference):
            raise ReferenceError("ReferenceBundle.annotation must be AnnotationReference")
        if not isinstance(self.coordinate_contract, CoordinateContract):
            raise ReferenceError(
                "ReferenceBundle.coordinate_contract must be CoordinateContract"
            )
        if not isinstance(self.sequence_extraction_ready, bool):
            raise ReferenceError("sequence_extraction_ready must be boolean")
        if not isinstance(self.sequence_extraction_blockers, (tuple, list)):
            raise ReferenceError("sequence_extraction_blockers must be an array")
        sequence_blockers = tuple(
            _as_string(
                blocker,
                f"ReferenceBundle.sequence_extraction_blockers[{index}]",
            )
            for index, blocker in enumerate(self.sequence_extraction_blockers)
        )
        if self.sequence_extraction_ready and sequence_blockers:
            raise ReferenceError(
                "a sequence-ready reference may not retain extraction blockers"
            )
        if not self.sequence_extraction_ready and not sequence_blockers:
            raise ReferenceError(
                "a sequence-unready reference must declare extraction blockers"
            )
        _set_frozen(self, "sequence_extraction_blockers", sequence_blockers)
        if not isinstance(self.prohibited_paths, (tuple, list)):
            raise ReferenceError("ReferenceBundle.prohibited_paths must be an array")
        prohibited_paths = tuple(
            _as_string(path, f"ReferenceBundle.prohibited_paths[{index}]")
            for index, path in enumerate(self.prohibited_paths)
        )
        if len(set(prohibited_paths)) != len(prohibited_paths):
            raise ReferenceError("ReferenceBundle.prohibited_paths contains duplicates")
        _set_frozen(self, "prohibited_paths", prohibited_paths)
        _set_frozen(
            self,
            "prohibited_reason",
            _as_string(self.prohibited_reason, "ReferenceBundle.prohibited_reason"),
        )
        missing_known = KNOWN_TRUNCATED_GENCODE_V49_PATHS - set(prohibited_paths)
        if missing_known:
            raise ReferenceError(
                "prohibited reference registry omits known truncated GENCODE v49 path(s): "
                + ", ".join(sorted(missing_known))
            )
        if self._base_dir is not None:
            _set_frozen(self, "_base_dir", Path(self._base_dir).resolve())
        if self.annotation_release.casefold().replace(" ", "") in {
            "v49",
            "gencodev49",
            "gencode_v49",
        } and self.annotation.sha256 != KNOWN_GENCODE_V49_ANNOTATION_SHA256:
            raise ReferenceError(
                "GENCODE v49 annotation must use the project-verified full-file SHA-256"
            )
        if (
            self.assembly == "GRCh38"
            and self.assembly_patch == "p14"
            and self.genome.sha256 != KNOWN_GRCH38_P14_GENOME_SHA256
        ):
            raise ReferenceError(
                "GRCh38.p14 genome must use the project-verified full-file SHA-256"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ReferenceBundle":
        root = _strict_mapping(
            value,
            allowed=cls._ROOT_FIELDS,
            required=cls._ROOT_REQUIRED,
            label="ReferenceBundle document",
        )
        reference = _strict_mapping(
            root["reference"],
            allowed=cls._REFERENCE_FIELDS,
            required=cls._REFERENCE_FIELDS,
            label="ReferenceBundle.reference",
        )
        prohibited = _strict_mapping(
            root["prohibited_references"],
            allowed=cls._PROHIBITED_FIELDS,
            required=cls._PROHIBITED_FIELDS,
            label="ReferenceBundle.prohibited_references",
        )
        return cls(
            bundle_id=reference["bundle_id"],
            assembly=reference["assembly"],
            assembly_patch=reference["assembly_patch"],
            annotation_release=reference["annotation_release"],
            genome=GenomeReference.from_dict(reference["genome"]),
            annotation=AnnotationReference.from_dict(reference["annotation"]),
            sequence_extraction_ready=reference["sequence_extraction_ready"],
            sequence_extraction_blockers=tuple(
                reference["sequence_extraction_blockers"]
            ),
            coordinate_contract=CoordinateContract.from_dict(
                root["coordinate_contract"]
            ),
            prohibited_paths=tuple(prohibited["paths"]),
            prohibited_reason=prohibited["reason"],
        )

    @classmethod
    def load_toml(cls, path: str | Path) -> "ReferenceBundle":
        source = Path(path)
        try:
            document = source.read_text(encoding="utf-8")
        except OSError as error:
            raise ReferenceError(f"cannot read {source}: {error}") from error
        bundle = cls.from_toml(document)
        _set_frozen(bundle, "_base_dir", Path(path).resolve().parent)
        return bundle

    def to_dict(self) -> dict[str, Any]:
        return {
            "reference": {
                "bundle_id": self.bundle_id,
                "assembly": self.assembly,
                "assembly_patch": self.assembly_patch,
                "annotation_release": self.annotation_release,
                "genome": self.genome.to_dict(),
                "annotation": self.annotation.to_dict(),
                "sequence_extraction_ready": self.sequence_extraction_ready,
                "sequence_extraction_blockers": list(
                    self.sequence_extraction_blockers
                ),
            },
            "coordinate_contract": self.coordinate_contract.to_dict(),
            "prohibited_references": {
                "paths": list(self.prohibited_paths),
                "reason": self.prohibited_reason,
            },
        }

    def validate(
        self,
        root: str | Path | None = None,
        *,
        require_indexed_sequence: bool = False,
    ) -> dict[str, str]:
        validation_root = root if root is not None else self._base_dir
        genome_path = self.genome.validate(validation_root)
        annotation_path = self.annotation.validate(validation_root)
        if (
            self.annotation_release.casefold().replace(" ", "")
            in {"v49", "gencodev49", "gencode_v49"}
            and annotation_path.name
            != "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
        ):
            raise ReferenceError(
                "GENCODE v49 annotation path is not the canonical full annotation filename"
            )
        genome_contigs = _fasta_contigs(genome_path)
        annotation_contigs = _annotation_contigs(annotation_path)
        missing_analysis_contigs = sorted(
            set(self.coordinate_contract.analysis_contigs) - genome_contigs
        )
        if missing_analysis_contigs:
            raise ReferenceError(
                "analysis coordinate contract contains contigs absent from the genome FASTA: "
                + ", ".join(missing_analysis_contigs)
            )
        missing_contigs = sorted(annotation_contigs - genome_contigs)
        if missing_contigs:
            preview = ", ".join(missing_contigs[:10])
            raise ReferenceError(
                "annotation contains contigs absent from the genome FASTA: " + preview
            )
        if require_indexed_sequence and not self.sequence_extraction_ready:
            raise ReferenceError(
                "indexed sequence extraction is admission-blocking: "
                + "; ".join(self.sequence_extraction_blockers)
            )
        return {
            "genome": genome_path.as_posix(),
            "annotation": annotation_path.as_posix(),
        }


__all__ = [
    "AnnotationReference",
    "CoordinateContract",
    "GenomeReference",
    "KNOWN_GENCODE_V49_ANNOTATION_SHA256",
    "KNOWN_GRCH38_P14_GENOME_SHA256",
    "KNOWN_TRUNCATED_GENCODE_V49_PATHS",
    "ReferenceBundle",
    "ReferenceError",
]
