"""
The exposure data model.

A benchmark that ranks encoders without resolving pretraining exposure publishes
a ranking it cannot read. In this substrate the two highest-scoring encoders,
Geneformer and TranscriptFormer, have ZERO clean held-out studies: every study
is either inside their pretraining corpus or of unresolved provenance. Their
numbers are not comparable to anything.

The rule this module enforces is that such a model reports NO NUMBER AT ALL.
Not a number with an asterisk. Not an empty string in a numeric column. Not a
null that a downstream read turns into NaN and someone sorts by.

That is enforced here by TYPE, not by documentation:

  * a comparable result can only be constructed when at least one held-out
    study is clean -- the constructor raises otherwise;
  * a non-comparable result has NO numeric attribute, and reaching for one
    raises an error that names the model and says why;
  * the leaderboard refuses any entry that is not a comparable result;
  * there is NO function in this package that returns a bare score.

Exposure is resolved PER STUDY, never per model. The two disagree in real data:
geneformer_v2_316m declares checkpoint state `target_label_unexposed` while
being `encoder_seen` on three of the seven held-out studies.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

# Read from scripts/score_encoder_field_common_head_50000.py rather than
# reinvented. An unresolved corpus is NOT evidence of no exposure, so `unknown`
# is its own disposition and never folded into clean.
CLEAN_STATES = frozenset({"clean_declared", "target_label_unexposed"})
CONFOUNDED_STATES = frozenset({"encoder_seen", "continual_seen", "reference_only"})

#: Attribute names that would carry a comparable quantity. Reaching for any of
#: these on a NotComparable is the mistake this module exists to prevent.
NUMERIC_ATTRIBUTES = frozenset({
    "delta", "lower", "upper", "median", "estimate", "score", "value",
    "macro_f1", "point", "oriented_delta", "probability_model_better",
})


class NotComparableError(RuntimeError):
    """Raised when a number is requested from a model that has none."""


class NoCleanStudyError(ValueError):
    """Raised when a comparable result is constructed with no clean study."""


class ExposureLedgerError(ValueError):
    """Raised when the three dispositions do not partition the study roster."""


def disposition(state: str) -> str:
    """clean / confounded / unresolved for one per-study exposure state."""
    if state in CLEAN_STATES:
        return "clean"
    if state in CONFOUNDED_STATES:
        return "confounded"
    return "unresolved"


@dataclass(frozen=True)
class ExposureLedger:
    """A partition of the study roster for one model.

    The invariant is checked, not assumed: the three sets are pairwise disjoint
    and their union is exactly the roster. A study cannot be silently dropped,
    which is how a confounded study becomes invisible.
    """

    model_id: str
    clean: tuple[str, ...]
    confounded: tuple[str, ...]
    unresolved: tuple[str, ...]
    roster: tuple[str, ...]

    def __post_init__(self) -> None:
        c, f, u = set(self.clean), set(self.confounded), set(self.unresolved)
        roster = set(self.roster)
        union = c | f | u
        if union != roster:
            raise ExposureLedgerError(
                f"{self.model_id}: the three dispositions do not cover the roster. "
                f"missing from the partition: {sorted(roster - union)}; "
                f"not in the roster: {sorted(union - roster)}"
            )
        overlaps = (c & f) | (c & u) | (f & u)
        if overlaps:
            raise ExposureLedgerError(
                f"{self.model_id}: studies in more than one disposition: {sorted(overlaps)}"
            )

    @classmethod
    def from_per_study(cls, model_id: str, per_study: Mapping[str, str]) -> "ExposureLedger":
        buckets: dict[str, list[str]] = {"clean": [], "confounded": [], "unresolved": []}
        for study in sorted(per_study):
            buckets[disposition(per_study[study])].append(study)
        return cls(
            model_id=model_id,
            clean=tuple(buckets["clean"]),
            confounded=tuple(buckets["confounded"]),
            unresolved=tuple(buckets["unresolved"]),
            roster=tuple(sorted(per_study)),
        )

    @property
    def has_clean_study(self) -> bool:
        return len(self.clean) > 0


@dataclass(frozen=True)
class ComparableDelta:
    """A paired delta against a named reference, on clean held-out studies only.

    Cannot exist without a clean study: the constructor raises. This is the only
    object in the package that carries a number.
    """

    model_id: str
    reference_id: str
    head_id: str
    metric: str
    clean_studies: tuple[str, ...]
    confounded_studies: tuple[str, ...]
    unresolved_studies: tuple[str, ...]
    delta: float
    lower: float
    upper: float
    median: float
    probability_model_better: float
    n_resamples: int
    resampling_unit: str = "donor"
    comparable: bool = field(default=True, init=False)

    def __post_init__(self) -> None:
        if not self.clean_studies:
            raise NoCleanStudyError(
                f"{self.model_id}: a comparable delta cannot be built with zero clean "
                f"held-out studies. Build a NotComparable instead."
            )

    @property
    def interval_crosses_zero(self) -> bool:
        return self.lower <= 0.0 <= self.upper

    @property
    def verdict(self) -> str:
        if self.interval_crosses_zero:
            return "interval_crosses_zero"
        return "model_above_reference" if self.delta > 0 else "model_below_reference"

    def to_record(self) -> dict:
        return {
            "model_id": self.model_id,
            "reference_id": self.reference_id,
            "head_id": self.head_id,
            "metric": self.metric,
            "comparable": True,
            "clean_studies": ",".join(self.clean_studies),
            "clean_study_count": len(self.clean_studies),
            "confounded_studies": ",".join(self.confounded_studies),
            "unresolved_studies": ",".join(self.unresolved_studies),
            "resampling_unit": self.resampling_unit,
            "n_resamples": self.n_resamples,
            "delta": self.delta,
            "lower": self.lower,
            "median": self.median,
            "upper": self.upper,
            "probability_model_better": self.probability_model_better,
            "interval_crosses_zero": self.interval_crosses_zero,
            "verdict": self.verdict,
        }


@dataclass(frozen=True)
class NotComparable:
    """A model with no clean held-out study.

    Carries no number and cannot be made to produce one. Reaching for `.delta`,
    `.lower`, `.estimate` or any other numeric name raises NotComparableError
    with the model named and the reason given, so the failure teaches instead of
    quietly returning None.
    """

    model_id: str
    reference_id: str
    head_id: str
    metric: str
    reason: str
    confounded_studies: tuple[str, ...]
    unresolved_studies: tuple[str, ...]
    comparable: bool = field(default=False, init=False)

    def __getattr__(self, name: str):
        # Reached only when normal attribute lookup has already failed, so the
        # declared fields above are unaffected.
        if name in NUMERIC_ATTRIBUTES:
            raise NotComparableError(
                f"{self.model_id} has no comparable {name}: {self.reason}. "
                f"confounded on [{', '.join(self.confounded_studies) or 'none'}]; "
                f"unresolved on [{', '.join(self.unresolved_studies) or 'none'}]. "
                f"There is no clean held-out study, so no number exists to report."
            )
        raise AttributeError(name)

    def to_record(self) -> dict:
        # The numeric keys are ABSENT, not null. A consumer that reads this
        # record cannot find a field to sort by.
        return {
            "model_id": self.model_id,
            "reference_id": self.reference_id,
            "head_id": self.head_id,
            "metric": self.metric,
            "comparable": False,
            "reason": self.reason,
            "clean_study_count": 0,
            "confounded_studies": ",".join(self.confounded_studies),
            "unresolved_studies": ",".join(self.unresolved_studies),
        }


Result = ComparableDelta | NotComparable


def build_result(
    *,
    model_id: str,
    reference_id: str,
    head_id: str,
    metric: str,
    ledger: ExposureLedger,
    delta: float | None = None,
    lower: float | None = None,
    upper: float | None = None,
    median: float | None = None,
    probability_model_better: float | None = None,
    n_resamples: int = 0,
) -> Result:
    """The only way to make a result. Exposure decides which type comes back."""
    if not ledger.has_clean_study:
        return NotComparable(
            model_id=model_id,
            reference_id=reference_id,
            head_id=head_id,
            metric=metric,
            reason="no_clean_held_out_study_exists_for_this_encoder",
            confounded_studies=ledger.confounded,
            unresolved_studies=ledger.unresolved,
        )
    if delta is None:
        raise ValueError(f"{model_id}: clean studies exist but no delta was computed")
    return ComparableDelta(
        model_id=model_id,
        reference_id=reference_id,
        head_id=head_id,
        metric=metric,
        clean_studies=ledger.clean,
        confounded_studies=ledger.confounded,
        unresolved_studies=ledger.unresolved,
        delta=float(delta),
        lower=float(lower),
        upper=float(upper),
        median=float(median),
        probability_model_better=float(probability_model_better),
        n_resamples=int(n_resamples),
    )


def partition(results: Iterable[Result]) -> tuple[list[ComparableDelta], list[NotComparable]]:
    comparable, withheld = [], []
    for r in results:
        (withheld if isinstance(r, NotComparable) else comparable).append(r)
    return comparable, withheld


def leaderboard(entries: Sequence[Result]) -> list[ComparableDelta]:
    """Rank comparable results. Refuses anything else.

    Passing a NotComparable here is the mistake that would publish an
    uninterpretable ranking, so it raises rather than being skipped.
    """
    for e in entries:
        if isinstance(e, NotComparable):
            raise NotComparableError(
                f"refusing to rank {e.model_id}: {e.reason}. A model with no clean "
                f"held-out study has no place in a leaderboard. Use partition() to "
                f"separate it out and report it in its own table."
            )
        if not isinstance(e, ComparableDelta):
            raise TypeError(f"leaderboard() takes ComparableDelta, got {type(e).__name__}")
    return sorted(entries, key=lambda e: e.delta, reverse=True)


LEADERBOARD_FIELDS = [
    "model_id", "reference_id", "head_id", "metric", "clean_studies",
    "clean_study_count", "confounded_studies", "unresolved_studies",
    "resampling_unit", "n_resamples", "delta", "lower", "median", "upper",
    "probability_model_better", "interval_crosses_zero", "verdict",
]

NOT_COMPARABLE_FIELDS = [
    "model_id", "reference_id", "head_id", "metric", "reason",
    "clean_study_count", "confounded_studies", "unresolved_studies",
]


def _write_tsv(path, fields: Sequence[str], records: Sequence[Mapping]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(fields) + "\n")
        for rec in records:
            fh.write("\t".join(str(rec[f]) for f in fields) + "\n")


def write_leaderboard(path, entries: Sequence[Result]) -> list[ComparableDelta]:
    """Write the ranking. Inherits leaderboard()'s refusal, so this file
    structurally cannot contain a model with zero clean held-out studies."""
    ranked = leaderboard(entries)
    _write_tsv(path, LEADERBOARD_FIELDS, [e.to_record() for e in ranked])
    return ranked


def write_not_comparable(path, entries: Sequence[NotComparable]) -> None:
    """Write the withheld models. No numeric column exists in this file."""
    bad = [f for f in NOT_COMPARABLE_FIELDS if f in NUMERIC_ATTRIBUTES]
    if bad:
        raise AssertionError(f"not_comparable.tsv must carry no numeric column, found {bad}")
    _write_tsv(path, NOT_COMPARABLE_FIELDS, [e.to_record() for e in entries])


def dump_results(path, results: Sequence[Result]) -> None:
    comparable, withheld = partition(results)
    payload = {
        "reportable_quantity": (
            "paired delta against the task's own native baseline, restricted to the "
            "model's clean held-out studies, with a donor cluster bootstrap interval. "
            "Never a bare score."
        ),
        "n_comparable": len(comparable),
        "n_withheld_for_exposure": len(withheld),
        "leaderboard": [e.to_record() for e in leaderboard(comparable)],
        "not_comparable": [e.to_record() for e in withheld],
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
