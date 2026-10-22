"""Restore legacy chromosome.position sentinels from native member coordinates.

The decimal encoding can lose trailing zeros. Neither fixed-width padding nor
the member mean identifies the original sentinel. Require one decimal extension
within the declared window of every member of the same credible set instead.
"""
from __future__ import annotations

import operator
import re

import numpy as np
import pandas as pd


def _chromosome(value: object) -> str:
    text = str(value).strip()
    if text.startswith("chr"):
        text = text[3:]
    if text in {"X", "Y"}:
        return text
    if not re.fullmatch(r"[1-9][0-9]*", text) or not 1 <= int(text) <= 22:
        raise ValueError(f"invalid chromosome: {value!r}")
    return str(int(text))


def _locus(value: object) -> tuple[str, int, str]:
    text = str(value).strip()
    match = re.fullmatch(r"(chr)?([1-9][0-9]*|X|Y)\.([0-9]+)", text)
    if match is None:
        raise ValueError(f"invalid chromosome.position locus: {value!r}")
    chrom = _chromosome(match.group(2))
    digits = match.group(3)
    encoded_pos = int(digits)
    if encoded_pos <= 0:
        raise ValueError(f"nonpositive encoded locus position: {value!r}")
    # Equivalent literal encodings may retain different numbers of trailing zeros.
    return chrom, encoded_pos, f"{chrom}.{digits.rstrip('0')}"


def recover_locus_positions(frame: pd.DataFrame, members: pd.DataFrame,
                            max_distance_bp: int = 1_000_000) -> pd.DataFrame:
    """Recover one positive integer sentinel per credible-set UID, or fail.

    ``members`` must include the full native member set, including coding
    members when recovering a noncoding-context subset. Metadata and coordinates
    are checked before recovery. The input frame, row order and index are kept.
    ``max_distance_bp`` is the original sentinel/member window, not a clustering
    sensitivity parameter; changing a clustering window must not change recovery.
    """
    try:
        distance = operator.index(max_distance_bp)
    except TypeError as exc:
        raise ValueError("max_distance_bp must be a nonnegative integer") from exc
    if isinstance(max_distance_bp, bool) or distance < 0:
        raise ValueError("max_distance_bp must be a nonnegative integer")
    for table, required, label in [
        (frame, {"credible_set_uid", "locus"}, "architecture/context"),
        (members, {"credible_set_uid", "chromosome", "position"}, "members"),
    ]:
        missing = required - set(table.columns)
        if missing:
            raise ValueError(f"{label} missing columns: {sorted(missing)}")
        if table.empty:
            raise ValueError(f"{label} is empty")
        if table["credible_set_uid"].isna().any() or (
                table["credible_set_uid"].astype(str).str.strip() == "").any():
            raise ValueError(f"{label} has missing credible-set identifiers")
    if frame["credible_set_uid"].duplicated().any():
        raise ValueError("credible_set_uid is not unique in the input frame")

    m = members.copy()
    m["credible_set_uid"] = m["credible_set_uid"].astype(str)
    m["native_chrom"] = m["chromosome"].map(_chromosome)
    positions = pd.to_numeric(m["position"], errors="coerce")
    valid = np.isfinite(positions) & (positions > 0) & (positions == np.floor(positions))
    if not valid.all() or (positions > np.iinfo(np.int64).max).any():
        raise ValueError("native member positions must be finite positive integers")
    m["native_pos"] = positions.astype("int64")
    if "locus" in m:
        m["native_locus"] = m["locus"].map(lambda value: _locus(value)[2])

    bounds = {}
    for uid, block in m.groupby("credible_set_uid", sort=False):
        chroms = block["native_chrom"].unique()
        if len(chroms) != 1:
            raise ValueError(f"mixed native chromosomes for credible set {uid}")
        loci = block["native_locus"].unique() if "native_locus" in block else []
        if len(loci) > 1:
            raise ValueError(f"inconsistent member locus metadata for {uid}")
        bounds[uid] = (chroms[0], int(block["native_pos"].min()),
                       int(block["native_pos"].max()), loci[0] if len(loci) else None)

    recovered = []
    for row in frame[["credible_set_uid", "locus"]].itertuples(index=False):
        uid = str(row.credible_set_uid)
        if uid not in bounds:
            raise ValueError(f"missing native members for credible set {uid}")
        chrom, encoded, locus_key = _locus(row.locus)
        native_chrom, low, high, member_locus = bounds[uid]
        if chrom != native_chrom:
            raise ValueError(f"locus/member chromosome mismatch for {uid}")
        if member_locus is not None and locus_key != member_locus:
            raise ValueError(f"locus/member metadata mismatch for {uid}")
        # Intersection of all member-centred windows, not distance to a mean.
        lower, upper = max(1, high - distance), low + distance
        candidate, zeros, supported = encoded, 0, []
        while candidate <= upper:
            if candidate >= lower:
                supported.append((candidate, zeros))
            candidate *= 10
            zeros += 1
        if len(supported) != 1:
            raise ValueError(f"credible set {uid}: expected one supported trailing-zero "
                             f"restoration, found {len(supported)}")
        pos, zeros = supported[0]
        recovered.append((chrom, pos, encoded, zeros, low, high))

    out = frame.copy()
    for i, column in enumerate(["chrom", "pos", "encoded_pos", "trailing_zeros_restored",
                                "member_min_position", "member_max_position"]):
        out[column] = [row[i] for row in recovered]
    return out
