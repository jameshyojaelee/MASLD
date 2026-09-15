"""Explicit multimodal pairing, missingness, and biological-unit requirements."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Any, Hashable, Iterable, Mapping


class TopologyError(ValueError):
    """Raised when assay topology or inferential units are misrepresented."""


PAIRING_LEVELS = frozenset(
    {
        "same_molecule",
        "same_cell",
        "same_nucleus",
        "same_section",
        "adjacent_section",
        "same_sample_different_aliquot",
        "same_donor_different_tissue",
        "same_study_unpaired",
    }
)
MISSING_STATES = frozenset(
    {
        "observed",
        "structurally_missing",
        "not_applicable",
        "below_qc",
        "unavailable_permission",
        "join_unresolved",
        "withheld_sealed",
        "derivable_not_processed",
    }
)
TECHNICAL_UNITS = frozenset({"cell", "spot", "array", "section", "run", "guide"})
BIOLOGICAL_UNITS = frozenset(
    {
        "donor",
        "participant",
        "animal",
        "sample_independent",
        "experimental_replicate",
        "independent_experimental_batch",
    }
)


@dataclass(frozen=True, slots=True)
class MaskedModality:
    state: str
    observed_mask: bool
    value: Any = None

    def validate(self) -> None:
        if self.state not in MISSING_STATES:
            raise TopologyError(f"unknown missingness state: {self.state}")
        if self.state == "observed":
            if not self.observed_mask or self.value is None:
                raise TopologyError("observed modalities require a true mask and a value")
        elif self.observed_mask or self.value is not None:
            raise TopologyError("missing modalities require a false mask and null value, never zero")


def assert_pairing(dataset_id: str, pairing_level: str) -> None:
    if pairing_level not in PAIRING_LEVELS:
        raise TopologyError(f"unknown pairing level: {pairing_level}")
    if dataset_id.upper() == "GSE244832" and pairing_level in {
        "same_molecule",
        "same_cell",
        "same_nucleus",
    }:
        raise TopologyError("GSE244832 supports same-donor pseudobulk, not cell-level pairing")


def assert_biological_replication(rows: Iterable[Mapping[str, Any]]) -> None:
    """Reject analyses whose declared replicate is a technical unit."""

    for row in rows:
        unit = str(row.get("unit_of_replication", ""))
        if unit in TECHNICAL_UNITS:
            raise TopologyError(f"{unit}s are not independent biological replicates")
        if unit not in BIOLOGICAL_UNITS:
            raise TopologyError(f"unrecognized biological replication unit: {unit}")


def donor_pseudobulk(
    rows: Iterable[Mapping[str, Any]],
    *,
    donor_key: str,
    feature_key: str,
    count_key: str,
) -> dict[Hashable, dict[Hashable, float]]:
    """Sum assay-native counts to donor level; no cell is treated as a replicate."""

    result: dict[Hashable, dict[Hashable, float]] = {}
    for row in rows:
        donor = row[donor_key]
        feature = row[feature_key]
        raw_value = row[count_key]
        if isinstance(raw_value, bool):
            raise TopologyError("pseudobulk count inputs must be numeric, not boolean")
        try:
            value = float(raw_value)
        except (TypeError, ValueError) as error:
            raise TopologyError("pseudobulk count inputs must be numeric") from error
        if not isfinite(value) or value < 0:
            raise TopologyError("pseudobulk count inputs must be finite and nonnegative")
        result.setdefault(donor, {})[feature] = result.setdefault(donor, {}).get(feature, 0.0) + value
    return result
