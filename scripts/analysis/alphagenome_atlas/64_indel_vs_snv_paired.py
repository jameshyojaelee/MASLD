#!/usr/bin/env python3
"""Step 64: are the indels the Atlas defers predicted to do MORE than the SNVs it serves at the same signals?

Each indel is paired against the SNV controls drawn at its own signals, on absolute effects, with the 1-Mb
block as the resampling unit. The family is the five channels and BH is applied within each stratum.

Three strata, because 87.7% of this Resource's deferred indels are orientation-ambiguous: a VCF indel's short
allele is a prefix of its long one, so whenever the long allele is at the reference both readings are
representable and the reference cannot say which. For those, the alternate sequence that was scored rests on
the source's allele order rather than on evidence. `dbsnp_oriented` is the subset whose canonical REF/ALT dbSNP supplies -- the primary stratum, since dbSNP
names the variant the GWAS actually tested. `dbsnp_unique_record` drops the rows step 68 could only orient
with its allele-frequency tie-break, so a result that depends on that weaker rule shows up as a difference
between the two. `reference_resolved` is the smaller subset the reference alone orients, kept as the most
conservative view; it is read from step 67, which is the authoritative reference-only answer.

Outputs (tables/): indel_vs_snv_paired.json
"""

from __future__ import annotations

import csv
import importlib.util
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la

HERE = pathlib.Path(__file__).resolve().parent
CHANNELS = ("rna", "atac", "dnase", "h3k27ac", "splice")
STRATA = ("all", "dbsnp_oriented", "dbsnp_unique_record", "reference_resolved")
N_BOOT = 2000
SEED = 20260913


def _rescue():
    spec = importlib.util.spec_from_file_location("indel_rescue", HERE / "63_indel_rescue.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def benjamini_hochberg(pvals) -> list:
    """BH q-values, order preserved; a NaN p stays NaN rather than becoming significant."""
    p = list(pvals)
    idx = [i for i, v in enumerate(p) if v == v]
    m = len(idx)
    q = [float("nan")] * len(p)
    if not m:
        return q
    order = sorted(idx, key=lambda i: p[i])
    running = 1.0
    for rank, i in reversed(list(enumerate(order, start=1))):
        running = min(running, float(p[i]) * m / rank)
        q[i] = running
    return q


def stratum(rows, name: str) -> list:
    """Scored indels in one orientation stratum."""
    if name not in STRATA:
        raise la.ContractError(f"unknown stratum {name!r}; expected one of {STRATA}")
    ind = [r for r in rows if r.get("arm") == "indel" and r.get("state") == "scored"]
    if name == "all":
        return ind
    if name == "dbsnp_oriented":
        return [r for r in ind if str(r.get("orientation_source", "")) == "dbsnp"]
    if name == "dbsnp_unique_record":
        return [r for r in ind if str(r.get("orientation_source", "")) == "dbsnp"
                and str(r.get("dbsnp_rule", "")) == "unique_record"]
    return [r for r in ind if str(r.get("orientation_ambiguous", "")) == "False"]


def assert_distinct_variants(rows) -> None:
    """Refuse two scored indel rows for one variant (same hg38 site, same unordered allele pair).

    Two rows for one variant are one observation. The 152-target run carried 8 such pairs, 3 of them scored
    as both reciprocal variants, because each source id was its own target.
    """
    seen = {}
    for r in rows:
        if r.get("arm") != "indel" or r.get("state") != "scored":
            continue
        key = (r["chrom"], int(r["pos_hg38"]), frozenset((str(r["ref"]).upper(), str(r["alt"]).upper())))
        if key in seen:
            raise la.ContractError(f"one variant scored twice: {seen[key]} and {r.get('source_variant_id')} at {key[:2]}")
        seen[key] = r.get("source_variant_id")


def reference_resolved_ids(repair_rows) -> set:
    """Source ids the REFERENCE alone orients, read from step 67 rather than from the effects table.

    Step 63 runs before this and its `orientation_ambiguous` column was written by a version that zeroed
    the flag for dbSNP rows; step 67's table is the authoritative reference-only answer and needs no API
    call to re-derive.
    """
    return {r["source_variant_id"] for r in repair_rows
            if r.get("repair_state") == "recovered" and r.get("orientation_ambiguous") == "False"}


def channel_stats(rescue, indels, controls, channel: str) -> dict:
    pr = rescue.paired_rows(indels, controls, channel)
    st = rescue.block_sign_test(pr)
    vals = np.asarray([p["value"] for p in pr], dtype=float)
    by = {}
    for p in pr:
        by.setdefault(p["analysis_block"], []).append(float(p["value"]))
    blocks = [np.asarray(v) for v in by.values()]
    rng = np.random.default_rng(SEED)
    boot = [float(np.median([float(blocks[i][rng.integers(blocks[i].size)])
                             for i in rng.integers(len(blocks), size=len(blocks))]))
            for _ in range(N_BOOT)] if blocks else []
    return {"n_paired_indels": len(pr), "n_blocks": st["n_blocks"],
            "median_paired_diff_log2": float(np.median(vals)) if vals.size else None,
            "lo": float(np.quantile(boot, 0.025)) if boot else None,
            "hi": float(np.quantile(boot, 0.975)) if boot else None,
            "blocks_positive": f"{st['n_blocks_positive']}/{st['n_blocks_nonzero']}",
            "sign_test_p": st["p_two_sided"]}


def main() -> None:
    rescue = _rescue()
    tables = la.out_root() / "tables"
    rows = list(csv.DictReader((tables / "indel_rescue_effects.tsv").open(), delimiter="\t"))
    assert_distinct_variants(rows)
    ctl = [r for r in rows if r.get("arm") == "snv_control"]
    out = {"n_controls": len(ctl), "family": "the five channels, within each stratum", "seed": SEED,
           "why_two_strata": ("a VCF indel's short allele is a prefix of its long one, so when the long allele "
                              "is at the reference both readings are representable; the scored alternate "
                              "sequence then rests on the source's allele order, not on the reference")}
    ref_ok = set()
    rp = sorted((la.out_root().parent).glob("p0-indel-uid-repair-*/tables/deferred_indel_uid_repair.tsv"))
    rp = [x for x in rp if not (x.parent.parent / "SUPERSEDED.txt").exists()]
    if rp:
        ref_ok = reference_resolved_ids(la.read_tsv(rp[-1]))
    for name in STRATA:
        ind = stratum(rows, name)
        if name == "reference_resolved" and ref_ok:
            ind = [r for r in stratum(rows, "all") if r.get("source_variant_id") in ref_ok]
        stats = {ch: channel_stats(rescue, ind, ctl, ch) for ch in CHANNELS}
        qs = benjamini_hochberg([stats[ch]["sign_test_p"] for ch in CHANNELS])
        for ch, q in zip(CHANNELS, qs):
            stats[ch]["bh_q"] = q
        out[name] = {"n_indels_scored": len(ind), **stats}
    json.dump(out, (tables / "indel_vs_snv_paired.json").open("w"), indent=1, default=float)
    for name in STRATA:
        la.log(f"step 64 [{name}]: n={out[name]['n_indels_scored']} " +
               " ".join(f"{ch} q={out[name][ch]['bh_q']:.2g}" for ch in CHANNELS))


if __name__ == "__main__":
    main()
