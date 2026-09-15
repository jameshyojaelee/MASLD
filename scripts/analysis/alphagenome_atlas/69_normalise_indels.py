#!/usr/bin/env python3
"""Step 69: canonical indel normalisation, so dbSNP's repeat representation can be matched.

The defect this closes. dbSNP collapses a homopolymer or short tandem repeat into ONE multi-allelic record
whose REF spans the whole repeat and whose ALTs enumerate the other run lengths: `C/CTT` at chr12:93522321
sits inside rs3038218 as REF=CTTTTTTT with twenty ALTs. Step 68 matched (REF, ALT) as strings, so it missed
every such variant and recorded it as absent from dbSNP. A sample of 400 showed 77% of the 24,902 "absent"
indels do have their allele pair at the exact dbSNP position -- they were never absent, only differently
written.

`normalise_indel` reduces a representation to the canonical left-aligned, trimmed form (the vt/bcftools
rule): trim a shared suffix, left-extend through a repeat when the alleles still end alike, then trim a
shared prefix. Two representations describe the same event exactly when their normal forms agree.

This module is a library for step 68; it performs no I/O of its own.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la


def normalise_indel(chrom: str, pos: int, ref: str, alt: str, fetch) -> tuple:
    """Canonical (pos, ref, alt) for one representation of one variant.

    `fetch` takes zero-based half-open coordinates. The left-extension step is what makes two spellings of
    a repeat indel converge; without it `C+7T -> C+9T` and `C -> CTT` stay different strings forever.
    """
    r, a = str(ref).upper(), str(alt).upper()
    p = int(pos)
    if not r or not a:
        raise la.ContractError(f"empty allele at {chrom}:{pos} {ref}/{alt}")
    if r == a:
        raise la.ContractError(f"reference and alternate are identical at {chrom}:{pos}: {r}")
    while True:
        if len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
            r, a = r[:-1], a[:-1]
            continue
        if r[-1] == a[-1] and (len(r) == 1 or len(a) == 1) and p > 1:
            base = fetch(chrom, p - 2, p - 1).upper()
            if not base:
                break
            r, a, p = base + r[:-1], base + a[:-1], p - 1
            continue
        break
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r, a, p = r[1:], a[1:], p + 1
    return p, r, a


def orientation_from_normalised(records, chrom: str, pos: int, a1: str, a2: str, fetch) -> dict:
    """Which of our two alleles dbSNP calls the reference, comparing canonical forms rather than strings.

    `records` is an iterable of (pos, REF, [ALT, ...]). Each (REF, ALTi) is a separate variant and is
    normalised on its own. Both of our orientations normalised and both found means dbSNP holds the two
    reciprocal events at that repeat and cannot attribute the source id to one; that is reported.
    """
    ours = {}
    for ref, alt in ((a1, a2), (a2, a1)):
        try:
            ours[normalise_indel(chrom, pos, ref, alt, fetch)] = (str(ref).upper(), str(alt).upper())
        except la.ContractError:
            continue
    if not ours:
        return {"state": "absent", "ref": None, "alt": None}
    seen = set()
    for rpos, rref, ralts in records:
        for ralt in ralts:
            try:
                key = normalise_indel(chrom, rpos, rref, ralt, fetch)
            except la.ContractError:
                continue
            if key in ours:
                seen.add(ours[key])
    if len(seen) == 1:
        ref, alt = seen.pop()
        return {"state": "resolved", "ref": ref, "alt": alt}
    if len(seen) > 1:
        return {"state": "both_in_dbsnp", "ref": None, "alt": None}
    return {"state": "absent", "ref": None, "alt": None}


def frequency_by_orientation(records, chrom: str, pos: int, a1: str, a2: str, fetch, freqs) -> dict:
    """Max allele frequency of the dbSNP ALT matching each of our two candidate orientations.

    `records` is (pos, REF, [ALT,...]) and `freqs` the parallel per-ALT frequency list for each record, so a
    repeat record's twenty alleles each carry their own frequency and the two that match us can be compared.
    """
    ours = {}
    for ref, alt in ((a1, a2), (a2, a1)):
        try:
            ours[normalise_indel(chrom, pos, ref, alt, fetch)] = (str(ref).upper(), str(alt).upper())
        except la.ContractError:
            continue
    out = {v: None for v in ours.values()}
    for (rpos, rref, ralts), rfreq in zip(records, freqs):
        for i, ralt in enumerate(ralts):
            try:
                key = normalise_indel(chrom, rpos, rref, ralt, fetch)
            except la.ContractError:
                continue
            if key not in ours:
                continue
            f = rfreq[i] if i < len(rfreq) else None
            if f is not None:
                cur = out[ours[key]]
                out[ours[key]] = f if cur is None else max(cur, f)
    return out
