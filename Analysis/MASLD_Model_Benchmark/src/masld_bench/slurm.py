"""Validated Slurm rendering and resource accounting.

This module renders scripts; it does not submit them.  Submission is an
explicit CLI action that also requires the approved campaign-plan hash.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import re
import shlex
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


class SlurmError(ValueError):
    """Raised when a request violates the project Slurm contract."""


_PARTITION_LIMIT_SECONDS = {
    "cpu": 14 * 24 * 60 * 60,
    "bigmem": 7 * 24 * 60 * 60,
    "io": 3 * 24 * 60 * 60,
    "gpu": 30 * 24 * 60 * 60,
}
_PROFILE_FIELDS = frozenset(
    {
        "partition",
        "cpus",
        "memory_gb",
        "wall_time",
        "gpus",
        "accelerator",
        "qos",
        "admission_blocking",
        "blockers",
    }
)
_PROFILE_REQUIRED = frozenset({"partition", "cpus", "memory_gb", "wall_time"})


def parse_duration(value: str) -> int:
    """Parse Slurm D-HH:MM:SS or HH:MM:SS into seconds."""

    match = re.fullmatch(r"(?:(\d+)-)?(\d+):(\d{2}):(\d{2})", value)
    if not match:
        raise SlurmError(f"invalid Slurm duration: {value!r}")
    days, hours, minutes, seconds = (int(part or 0) for part in match.groups())
    if minutes >= 60 or seconds >= 60 or (days and hours >= 24):
        raise SlurmError(f"invalid Slurm duration fields: {value!r}")
    return (((days * 24) + hours) * 60 + minutes) * 60 + seconds


@dataclass(frozen=True, slots=True)
class ResourceProfile:
    profile_id: str
    partition: str
    cpus: int
    memory_gb: int
    wall_time: str
    gpus: int = 0
    accelerator: str | None = None
    qos: str | None = None
    admission_blocking: bool = False
    blockers: tuple[str, ...] = ()

    @classmethod
    def from_mapping(cls, profile_id: str, value: Mapping[str, Any]) -> "ResourceProfile":
        if not isinstance(profile_id, str) or not profile_id:
            raise SlurmError("resource profile identifier must be a non-empty string")
        if not isinstance(value, Mapping):
            raise SlurmError(f"resource profile {profile_id!r} must be a mapping")
        unknown = sorted(set(value).difference(_PROFILE_FIELDS))
        missing = sorted(_PROFILE_REQUIRED.difference(value))
        if unknown or missing:
            details = []
            if missing:
                details.append("missing=" + ",".join(missing))
            if unknown:
                details.append("unknown=" + ",".join(unknown))
            raise SlurmError(
                f"resource profile {profile_id!r} has invalid fields: "
                + "; ".join(details)
            )
        for field_name in ("cpus", "memory_gb", "gpus"):
            raw = value.get(field_name, 0)
            if isinstance(raw, bool) or not isinstance(raw, int):
                raise SlurmError(
                    f"resource profile {profile_id!r} {field_name} must be an integer"
                )
        for field_name in ("partition", "wall_time"):
            if not isinstance(value[field_name], str) or not value[field_name]:
                raise SlurmError(
                    f"resource profile {profile_id!r} {field_name} must be a non-empty string"
                )
        for field_name in ("accelerator", "qos"):
            if field_name in value and value[field_name] is not None and (
                not isinstance(value[field_name], str) or not value[field_name]
            ):
                raise SlurmError(
                    f"resource profile {profile_id!r} {field_name} must be null or a non-empty string"
                )
        if "admission_blocking" in value and not isinstance(
            value["admission_blocking"], bool
        ):
            raise SlurmError(
                f"resource profile {profile_id!r} admission_blocking must be boolean"
            )
        raw_blockers = value.get("blockers", [])
        if not isinstance(raw_blockers, list) or any(
            not isinstance(item, str) or not item.strip() for item in raw_blockers
        ):
            raise SlurmError(
                f"resource profile {profile_id!r} blockers must be non-empty strings"
            )
        return cls(
            profile_id=profile_id,
            partition=str(value["partition"]),
            cpus=value["cpus"],
            memory_gb=value["memory_gb"],
            wall_time=str(value["wall_time"]),
            gpus=value.get("gpus", 0),
            accelerator=(str(value["accelerator"]) if value.get("accelerator") else None),
            qos=(str(value["qos"]) if value.get("qos") else None),
            admission_blocking=bool(value.get("admission_blocking", False)),
            blockers=tuple(str(item) for item in value.get("blockers", [])),
        )

    @classmethod
    def from_serialized(
        cls, profile_id: str, value: Mapping[str, Any]
    ) -> "ResourceProfile":
        """Reconstruct one profile emitted by :meth:`as_dict` without schema drift."""

        if not isinstance(value, Mapping):
            raise SlurmError(f"serialized resource profile {profile_id!r} must be a mapping")
        expected = _PROFILE_FIELDS | {"profile_id"}
        if set(value) != expected:
            missing = sorted(expected.difference(value))
            unknown = sorted(set(value).difference(expected))
            details = []
            if missing:
                details.append("missing=" + ",".join(missing))
            if unknown:
                details.append("unknown=" + ",".join(unknown))
            raise SlurmError(
                f"serialized resource profile {profile_id!r} has invalid fields: "
                + "; ".join(details)
            )
        if value["profile_id"] != profile_id:
            raise SlurmError(
                f"serialized resource profile key {profile_id!r} differs from its profile_id"
            )
        profile = cls.from_mapping(
            profile_id,
            {key: item for key, item in value.items() if key != "profile_id"},
        )
        if profile.as_dict() != dict(value):
            raise SlurmError(
                f"serialized resource profile {profile_id!r} is not canonical"
            )
        return profile

    @property
    def wall_seconds(self) -> int:
        return parse_duration(self.wall_time)

    def validate(self, *, for_scheduling: bool = False) -> None:
        if not isinstance(self.profile_id, str) or not re.fullmatch(
            r"[A-Za-z0-9_.-]+", self.profile_id
        ):
            raise SlurmError(f"invalid resource profile identifier: {self.profile_id!r}")
        if self.partition not in _PARTITION_LIMIT_SECONDS:
            raise SlurmError(f"unsupported partition: {self.partition}")
        if (
            isinstance(self.cpus, bool)
            or not isinstance(self.cpus, int)
            or not 1 <= self.cpus <= 600
        ):
            raise SlurmError(f"invalid CPU count for {self.profile_id}: {self.cpus}")
        if (
            isinstance(self.memory_gb, bool)
            or not isinstance(self.memory_gb, int)
            or self.memory_gb <= 0
        ):
            raise SlurmError(f"memory must be positive for {self.profile_id}")
        if self.wall_seconds <= 0:
            raise SlurmError("Slurm wall time must be positive")
        if self.wall_seconds > _PARTITION_LIMIT_SECONDS[self.partition]:
            raise SlurmError(
                f"{self.wall_time} exceeds {self.partition} partition maximum"
            )
        if self.partition == "cpu" and self.memory_gb > 210:
            raise SlurmError("requests above 210 GB must use bigmem")
        if (
            isinstance(self.gpus, bool)
            or not isinstance(self.gpus, int)
            or self.gpus < 0
        ):
            raise SlurmError(f"invalid GPU count for {self.profile_id}: {self.gpus}")
        if self.partition != "gpu" and self.gpus:
            raise SlurmError("GPU requests must use the gpu partition")
        if self.partition == "gpu" and self.gpus < 1:
            raise SlurmError("gpu partition profiles must request at least one GPU")
        if self.gpus and not self.accelerator:
            raise SlurmError("GPU profiles must name the tested accelerator runtime")
        if self.qos is not None and not re.fullmatch(r"[A-Za-z0-9_.-]+", self.qos):
            raise SlurmError(f"invalid QOS identifier: {self.qos!r}")
        if self.qos == "interactive":
            raise SlurmError("interactive QOS requires explicit per-session user authorization")
        if not isinstance(self.admission_blocking, bool):
            raise SlurmError("admission_blocking must be boolean")
        if any(not isinstance(item, str) or not item.strip() for item in self.blockers):
            raise SlurmError("resource blockers must be non-empty strings")
        if self.admission_blocking and not self.blockers:
            raise SlurmError("an admission-blocking resource profile must state its blockers")
        if for_scheduling and self.admission_blocking:
            reason = "; ".join(self.blockers) or "profile has unresolved admission checks"
            raise SlurmError(f"resource profile {self.profile_id} is admission-blocking: {reason}")
        if for_scheduling and self.accelerator and self.accelerator.upper() == "UNRESOLVED":
            raise SlurmError(f"resource profile {self.profile_id} has no frozen accelerator")

    def as_dict(self) -> dict[str, Any]:
        self.validate()
        return {
            "profile_id": self.profile_id,
            "partition": self.partition,
            "cpus": self.cpus,
            "memory_gb": self.memory_gb,
            "wall_time": self.wall_time,
            "gpus": self.gpus,
            "accelerator": self.accelerator,
            "qos": self.qos,
            "admission_blocking": self.admission_blocking,
            "blockers": list(self.blockers),
        }


def render_sbatch(
    *,
    job_name: str,
    profile: ResourceProfile,
    command: Sequence[str],
    stdout_path: str | Path,
    stderr_path: str | Path,
    environment: Mapping[str, str] | None = None,
    modules: Sequence[str] = (),
    working_directory: str | Path | None = None,
) -> str:
    profile.validate(for_scheduling=True)
    if not command:
        raise SlurmError("an empty job command is not permitted")
    if any("--array" in part for part in command):
        raise SlurmError("array directives must be represented in the reviewed campaign plan")
    log_paths: dict[str, Path] = {}
    for label, raw_path in (("stdout_path", stdout_path), ("stderr_path", stderr_path)):
        path = Path(raw_path)
        rendered = path.as_posix()
        if not path.is_absolute() or any(
            character in rendered for character in ("\n", "\r", "\0")
        ):
            raise SlurmError(f"{label} must be an absolute control-character-free path")
        log_paths[label] = path
    safe_job_name = re.sub(r"[^A-Za-z0-9_.-]", "_", job_name)[:128]
    directives = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={safe_job_name}",
        f"#SBATCH --partition={profile.partition}",
        f"#SBATCH --cpus-per-task={profile.cpus}",
        f"#SBATCH --mem={profile.memory_gb}G",
        f"#SBATCH --time={profile.wall_time}",
        f"#SBATCH --output={log_paths['stdout_path'].as_posix()}",
        f"#SBATCH --error={log_paths['stderr_path'].as_posix()}",
    ]
    if profile.gpus:
        accelerator = str(profile.accelerator).casefold()
        if not re.fullmatch(r"[a-z0-9_.-]+", accelerator) or accelerator in {
            "none",
            "unresolved",
        }:
            raise SlurmError(
                f"invalid schedulable accelerator for {profile.profile_id}: "
                f"{profile.accelerator!r}"
            )
        directives.append(f"#SBATCH --gres=gpu:{accelerator}:{profile.gpus}")
    if profile.qos:
        directives.append(f"#SBATCH --qos={profile.qos}")
    body = ["", "set -euo pipefail", "umask 027"]
    if working_directory is not None:
        workdir = Path(working_directory)
        if not workdir.is_absolute() or "\n" in workdir.as_posix():
            raise SlurmError("working_directory must be an absolute path")
        body.append(f"cd -- {shlex.quote(workdir.as_posix())}")
    for module in modules:
        if not re.fullmatch(r"[A-Za-z0-9._/+:-]+", str(module)):
            raise SlurmError(f"invalid module identifier: {module!r}")
        body.append(f"module load {shlex.quote(str(module))}")
    for key, value in sorted((environment or {}).items()):
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
            raise SlurmError(f"invalid environment variable name: {key!r}")
        body.append(f"export {key}={shlex.quote(value)}")
    body.extend(["", shlex.join([str(part) for part in command]), ""])
    return "\n".join([*directives, *body])


def resource_totals(
    profiles: Mapping[str, ResourceProfile], run_profile_ids: Iterable[str]
) -> dict[str, Any]:
    job_count = 0
    cpu_hours = Decimal(0)
    gpu_hours = Decimal(0)
    memory_gb_hours = Decimal(0)
    total_cpus = 0
    total_gpus = 0
    total_memory_gb = 0
    by_profile: dict[str, int] = {}
    for profile_id in run_profile_ids:
        if profile_id not in profiles:
            raise SlurmError(f"unknown resource profile: {profile_id}")
        profile = profiles[profile_id]
        profile.validate(for_scheduling=True)
        hours = Decimal(profile.wall_seconds) / Decimal(3600)
        job_count += 1
        cpu_hours += Decimal(profile.cpus) * hours
        gpu_hours += Decimal(profile.gpus) * hours
        memory_gb_hours += Decimal(profile.memory_gb) * hours
        total_cpus += profile.cpus
        total_gpus += profile.gpus
        total_memory_gb += profile.memory_gb
        by_profile[profile_id] = by_profile.get(profile_id, 0) + 1
    return {
        "jobs": job_count,
        "cpu_hours_requested": float(cpu_hours),
        "gpu_hours_requested": float(gpu_hours),
        "memory_gb_hours_requested": float(memory_gb_hours),
        "cpus_if_all_concurrent": total_cpus,
        "gpus_if_all_concurrent": total_gpus,
        "memory_gb_if_all_concurrent": total_memory_gb,
        "jobs_by_profile": dict(sorted(by_profile.items())),
    }


def load_resource_profiles(document: Mapping[str, Any]) -> dict[str, ResourceProfile]:
    raw_profiles = document.get("profiles", document.get("profile", {}))
    if not isinstance(raw_profiles, Mapping):
        raise SlurmError("resources.toml must contain a [profiles.*] table")
    profiles = {
        str(profile_id): ResourceProfile.from_mapping(str(profile_id), value)
        for profile_id, value in raw_profiles.items()
    }
    if not profiles:
        raise SlurmError("no resource profiles are defined")
    for profile in profiles.values():
        profile.validate()
    return profiles
