#!/usr/bin/env python3
"""Step 63: score, through the model API, the indels the Atlas point query refuses.

P6b found that the Atlas annotation's real coverage gap is a variant CLASS, not an ancestry. Across the
989-signal portfolio, 207,254 variant-signal pairs carrying 49.5 units of posterior mass are deferred
indels; 13 signals (10 of them direct MASLD traits) carry essentially all of their posterior on indels, and
the deferred set includes **HSD17B13 rs72613567** (hg19 chr4:88231392 T>TA), the best-established
protective regulatory variant in MASLD.

Neither Atlas route reaches them: point queries return gRPC UNIMPLEMENTED, and interval saturation returns
only the three substitutions per position (0 multi-base alleles in 104,352 variants across 40 archives,
checked 2026-09-13). The model API is a different service and takes an arbitrary sequence, so the indel can
be built into the window directly. `predict_sequence` requires exactly 1 Mb, so the alternate window is
length-compensated at its far end, ~500 kb from the readout.

What this can and cannot say: the model API returns a predicted TRACK, not a calibrated quantile, so these
effects are NOT comparable to Atlas quantiles and are never pooled with them. Each indel is reported beside
the reference prediction it is a ratio to, and beside a matched set of SNVs at the same locus scored the
same way, which is the only thing that makes the magnitude readable.

Outputs (tables/): indel_rescue_effects.tsv, indel_rescue_summary.json
"""

from __future__ import annotations

import json
import math
import re
from collections import defaultdict

import numpy as np

import lib_atlas as la

SOURCE_RE = re.compile(r"^(?P<assembly>[^:]+):(?P<chrom>chr[0-9XYM]+|[0-9]{1,2}|[XYM]):(?P<pos>\d+):(?P<ref>[ACGTacgt]+):(?P<alt>[ACGTacgt]+)$")


def parse_source_id(source_id: str) -> dict:
    """`hg19:chr4:88231392:T:TA` and `hg19:4:88231392:T:TA` are both used in signal_variant_weights."""
    m = SOURCE_RE.match(source_id.strip())
    if not m:
        raise la.ContractError(f"unparseable source variant id: {source_id!r}")
    chrom = m.group("chrom")
    return {"assembly": m.group("assembly"), "chrom": chrom if chrom.startswith("chr") else f"chr{chrom}",
            "pos": int(m.group("pos")), "ref": m.group("ref").upper(), "alt": m.group("alt").upper()}


def indel_class(ref: str, alt: str) -> str:
    if len(alt) > len(ref):
        return "insertion"
    if len(alt) < len(ref):
        return "deletion"
    return "snv"


def length_delta(ref: str, alt: str) -> int:
    return len(alt) - len(ref)


def block_of(chrom: str, pos: int) -> str:
    return f"{chrom}:{int(pos) // 1_000_000}"


def _f(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def paired_rows(indels: list, controls: list, channel: str) -> list:
    """Per-indel |effect| minus the median |effect| of the SNV controls drawn at its OWN signals.

    A median ratio over all indels against all controls compares two differently-composed sets; this pairs
    within signal. Absolute effects, because an indel and an SNV need not push the same way and a signed
    comparison would cancel rather than measure magnitude.
    """
    by = defaultdict(list)
    for c in controls:
        v = abs(_f(c.get(f"{channel}_log2")))
        if v == v:
            by[str(c.get("control_for", ""))].append(v)
    out = []
    for r in indels:
        v = abs(_f(r.get(f"{channel}_log2")))
        ctl = by.get(str(r.get("source_variant_id", "")), [])
        if v != v or not ctl:
            continue
        out.append({"source_variant_id": r.get("source_variant_id"),
                    "analysis_block": block_of(r["chrom"], r["pos_hg38"]),
                    "value": float(v - float(np.median(ctl))), "n_controls": len(ctl)})
    return out


def block_sign_test(rows: list) -> dict:
    """Sign test with the 1-Mb block as the unit: paired indels at one locus are not independent."""
    from math import comb
    by = defaultdict(list)
    for r in rows:
        by[r["analysis_block"]].append(float(r["value"]))
    meds = [float(np.median(v)) for _, v in sorted(by.items())]
    n = sum(1 for m in meds if m != 0.0)
    k = sum(1 for m in meds if m > 0)
    p = (sum(comb(n, i) for i in range(0, min(k, n - k) + 1)) * 2 / (2 ** n)) if n else float("nan")
    return {"n_blocks": len(meds), "n_blocks_positive": k, "n_blocks_nonzero": n,
            "median_block_effect": float(np.median(meds)) if meds else float("nan"),
            "p_two_sided": float(min(1.0, p)) if n else float("nan")}


def select_targets(mass: dict, min_mass: float, only: str, rep_of: dict | None = None) -> list:
    """Deferred indels worth a model call, heaviest first; `only` names an explicit subset for a probe.

    `rep_of` maps each source id to its variant group's representative, so an explicit id selects its
    group once whichever member was named.
    """
    if only.strip():
        wanted = [(rep_of or {}).get(v.strip(), v.strip()) for v in only.split(",") if v.strip()]
        return [v for v in dict.fromkeys(wanted) if v in mass]
    return sorted((v for v, m in mass.items() if m >= min_mass), key=lambda v: (-mass[v], v))


def variant_groups(indel_mass: dict, indel_sigs: dict, xw: dict) -> dict:
    """Source ids that name the SAME indel, grouped: representative -> members, summed mass, union of signals.

    Studies write one biallelic indel with its alleles in either order and with or without 'chr', so
    `T:TCA` and `TCA:T` at one site are one variant. Scoring each id as its own target counted 8 variants
    twice in the 152-target run and scored 3 of them as both reciprocal variants. The key is the hg38 site
    plus the UNORDERED allele pair; an id with no hg38 coordinate stays alone, since nothing establishes it
    is the same variant. The representative is the heaviest member (ties: smallest id).
    """
    keyed = defaultdict(list)
    for sid in indel_mass:
        c = xw.get(sid) or {}
        parsed = parse_source_id(sid)
        if c.get("hg38_chrom") and c.get("hg38_position_1based"):
            key = (c["hg38_chrom"], int(c["hg38_position_1based"]),
                   frozenset((parsed["ref"].upper(), parsed["alt"].upper())))
        else:
            key = ("", sid)
        keyed[key].append(sid)
    out = {}
    for members in keyed.values():
        members = sorted(members, key=lambda v: (-indel_mass[v], v))
        out[members[0]] = {"members": members, "mass": sum(indel_mass[v] for v in members),
                           "signals": set().union(*(indel_sigs[v] for v in members))}
    return out


def group_dbsnp_orientation(members: list, dbsnp: dict):
    """One dbSNP answer for ids naming the same variant; members with different answers are refused.

    The rule is 'unique_record' if any member was settled by a single record, since the frequency tie-break
    then only agreed with it.
    """
    hits = [dbsnp[m] for m in members if m in dbsnp]
    answers = {(h["dbsnp_ref"], h["dbsnp_alt"]) for h in hits}
    if len(answers) > 1:
        raise la.ContractError(f"dbSNP orients the ids of one variant differently: {members}")
    if not hits:
        return None
    ref, alt = answers.pop()
    rule = "unique_record" if any(h["dbsnp_rule"] == "unique_record" for h in hits) else "frequency"
    return {"dbsnp_ref": ref, "dbsnp_alt": alt, "dbsnp_rule": rule}


def fold_mask(mask, n_rows: int):
    """Fold a base-resolution mask onto a track's own row count.

    CHIP_HISTONE is 128-bp binned; applying a base mask to it selects the wrong rows without erroring,
    which is exactly how it went wrong in P5A.
    """
    m = np.asarray(mask, dtype=bool)
    if m.size == n_rows:
        return m
    if n_rows <= 0 or m.size % n_rows:
        raise la.ContractError(f"mask of {m.size} does not fold onto {n_rows} rows")
    fold = m.size // n_rows
    return m.reshape(n_rows, fold).any(axis=1)


def h3k27ac_columns(track_names) -> list:
    """CHIP_HISTONE carries 17 marks; only H3K27ac is the enhancer-activity readout used elsewhere here."""
    return [("H3K27AC" in str(n).upper()) for n in track_names]


def oriented_alleles(chrom: str, pos1: int, a1: str, a2: str, fetch) -> tuple:
    """hg38 reference and alternate for a source variant whose allele order is not guaranteed.

    The source ids come from GWAS summary statistics, where the two alleles are effect/other, not ref/alt.
    Testing only the first one against hg38 refused 12 of this lane's 152 targets on the first run, all of
    them written the other way round, and two of those carried the largest posterior mass in the list.
    `la.resolve_indel_alleles` tries the longer allele first, which is what makes the test discriminating
    for prefix-anchored indels.
    """
    r = la.resolve_indel_alleles(chrom, int(pos1), a1, a2, fetch)
    if r["variant_uid"] is None:
        raise la.ContractError(f"neither allele of {chrom}:{pos1} {a1}/{a2} is the hg38 reference")
    return r["hg38_ref"], r["hg38_alt"], r["allele_swap"], r["orientation_ambiguous"]


def orientation_for(source_id: str, chrom: str, pos1: int, a1: str, a2: str, fetch, dbsnp: dict) -> tuple:
    """Reference/alternate alleles for a target, dbSNP first, the reference rule second.

    dbSNP stores each variant in its canonical REF/ALT form and its rsIDs identify the variant the GWAS
    actually tested, so where dbSNP resolves an indel it outranks both the source's allele order and the
    reference-only rule. For a prefix-anchored indel that rule falls back on source order, which step 68b
    found right on 47.4% of 99,548 ambiguous indels with a study-convention truth, while dbSNP was right on
    99.86% (single record) and 99.85% (frequency tie-break). Where dbSNP is silent or holds both
    orientations at common frequency, the reference rule applies and the row is flagged so the caller can
    stratify on it.

    Fail-closed twice: a dbSNP record whose allele pair is not ours is refused rather than applied, and a
    dbSNP reference allele that is not at that position in our FASTA is refused rather than scored.
    """
    d = dbsnp.get(source_id)
    if d and d.get("dbsnp_ref") and d.get("dbsnp_alt"):
        ref, alt = d["dbsnp_ref"].upper(), d["dbsnp_alt"].upper()
        if {ref, alt} != {str(a1).upper(), str(a2).upper()}:
            raise la.ContractError(
                f"dbSNP record for {source_id} is {ref}/{alt}, not our pair {a1}/{a2}")
        obs = fetch(chrom, int(pos1) - 1, int(pos1) - 1 + len(ref)).upper()
        if obs != ref:
            raise la.ContractError(
                f"dbSNP reference {ref} for {source_id} is not at {chrom}:{pos1} (genome has {obs})")
        swap = "none" if ref == str(a1).upper() else "swapped"
        # The prefix ambiguity is a property of the alleles and the genome, not of dbSNP. Report it even
        # when dbSNP settles the orientation, or the conservative stratum collapses onto the dbSNP one and
        # stops being an independent check.
        ref_amb = la.resolve_indel_alleles(chrom, int(pos1), a1, a2, fetch)["orientation_ambiguous"]
        return ref, alt, swap, bool(ref_amb), "dbsnp"
    ref, alt, swap, amb = oriented_alleles(chrom, pos1, a1, a2, fetch)
    return ref, alt, swap, amb, "reference_or_source_order"


def dbsnp_orientations(rows) -> dict:
    """Source id -> dbSNP orientation, from step 68's RECONCILED answer only.

    Step 68 writes each assembly arm's call beside the reconciled `resolved_variant_uid`. Only the latter
    may be applied: an arm's own columns can hold a call the other arm refused. `dbsnp_rule` says whether a
    single dbSNP record decided it or the allele-frequency tie-break did, so the weaker rule can be
    stratified on.
    """
    out = {}
    for r in rows:
        uid = str(r.get("resolved_variant_uid") or "").strip()
        if not uid:
            continue
        parts = uid.split(":")
        if len(parts) != 4:
            raise la.ContractError(f"resolved_variant_uid is not chrom:pos:ref:alt: {uid!r}")
        rule = "frequency" if r.get("frequency_tie_break") == "resolved_by_frequency" else "unique_record"
        out[r["source_variant_id"]] = {"dbsnp_ref": parts[2], "dbsnp_alt": parts[3], "dbsnp_rule": rule}
    return out


def indel_alt_sequence(fetch, chrom: str, start0: int, end: int, pos1: int, ref: str, alt: str) -> str:
    """The alternate window, same length as the reference window, with the edit applied at pos1.

    An insertion pushes bases off the far end; a deletion pulls bases in from beyond it. The readout is
    local to the variant, which sits at the window centre, so the compensated edge is ~500 kb away.
    """
    width = end - start0
    i = pos1 - 1 - start0
    if i < 0 or i + len(ref) > width:
        raise la.ContractError(f"variant {chrom}:{pos1} {ref}>{alt} does not fit inside [{start0},{end})")
    extended = fetch(chrom, start0, end + max(0, len(ref) - len(alt))).upper()
    observed = extended[i:i + len(ref)]
    if observed != ref.upper():
        raise la.ContractError(f"reference mismatch at {chrom}:{pos1}: window has {observed}, variant says {ref}")
    out = extended[:i] + alt.upper() + extended[i + len(ref):]
    if len(out) < width:
        raise la.ContractError(f"cannot compensate window length for {chrom}:{pos1}")
    return out[:width]


def main() -> None:
    import csv

    import pysam
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    import atlas_query as aq

    ROOT = la.out_root()
    TABLES = ROOT / "tables"
    TRACK0 = la.track0_root() / "tables"
    FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
    WINDOW = dna_client.SEQUENCE_LENGTH_1MB
    OUTPUTS = ["RNA_SEQ", "ATAC", "DNASE", "CHIP_HISTONE", "SPLICE_SITE_USAGE"]
    LIVER_TERMS = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
    MIN_MASS = float(__import__("os").environ.get("AGA_INDEL_MIN_MASS", "0.01"))
    N_SNV_CONTROLS = int(__import__("os").environ.get("AGA_INDEL_CONTROLS", "6"))
    FLANK = 2_000
    EPS = 1e-9

    signals = {s["signal_uid"]: s for s in la.read_tsv(TRACK0 / "eligible_signals.tsv")}
    xw = {}
    with la.open_text(TRACK0 / "variant_crosswalk.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            xw[r["source_variant_id"]] = r

    # Deferred indels carrying real posterior mass, and the queried SNVs at the same signals as controls.
    indel_mass, indel_sigs, snv_pool = defaultdict(float), defaultdict(set), defaultdict(list)
    with la.open_text(TRACK0 / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            w = float(r["weight"])
            if r["exclusion_reason"] == "indel_deferred":
                indel_mass[r["source_variant_id"]] += w
                indel_sigs[r["source_variant_id"]].add(r["signal_uid"])
            elif r["mapping_status"] == "mapped" and r["in_query_set"] == "True" and w >= MIN_MASS:
                snv_pool[r["signal_uid"]].append((w, r["source_variant_id"]))

    groups = variant_groups(indel_mass, indel_sigs, xw)
    rep_of = {m: rep for rep, g in groups.items() for m in g["members"]}
    targets = select_targets({rep: g["mass"] for rep, g in groups.items()}, MIN_MASS,
                             __import__("os").environ.get("AGA_INDEL_ONLY", ""), rep_of)
    la.log(f"step 63: {len(indel_mass)} deferred indel ids naming {len(groups)} variants; "
           f"{len(targets)} variants carry >= {MIN_MASS} posterior mass")

    key = la.load_api_key()[0]
    model = dna_client.create(key)
    outs = [getattr(dna_client.OutputType, o) for o in OUTPUTS]
    fasta = pysam.FastaFile(FASTA)

    def summarise(output, mask_local) -> dict:
        out = {}
        for name, attr in (("rna", "rna_seq"), ("atac", "atac"), ("dnase", "dnase"),
                           ("h3k27ac", "chip_histone"), ("splice", "splice_site_usage")):
            td = getattr(output, attr, None)
            if td is None or td.values is None or td.values.size == 0:
                out[name] = math.nan
                continue
            v = np.asarray(td.values, dtype=float)
            keep = np.ones(v.shape[1], bool)
            if name == "h3k27ac":
                md = getattr(td, "metadata", None)
                if md is None or "name" not in md:
                    out[name] = math.nan
                    continue
                keep = np.asarray(h3k27ac_columns(md["name"].astype(str)), bool)
                if not keep.any():
                    out[name] = math.nan
                    continue
            m = fold_mask(mask_local, v.shape[0])
            out[name] = float(np.nanmean(np.abs(v[m, :][:, keep]))) if m.any() else math.nan
        return out

    def score(chrom: str, pos1: int, ref: str, alt: str, tag: str) -> dict:
        start0 = max(0, pos1 - 1 - WINDOW // 2)
        end = start0 + WINDOW
        ref_seq = fasta.fetch(chrom, start0, end).upper()
        if len(ref_seq) != WINDOW:
            return {}
        alt_seq = indel_alt_sequence(lambda c, a, b: fasta.fetch(c, a, b), chrom, start0, end, pos1, ref, alt)
        interval = genome.Interval(chromosome=chrom, start=start0, end=end)
        mask = np.zeros(WINDOW, bool)
        lo, hi = max(0, pos1 - 1 - start0 - FLANK), min(WINDOW, pos1 - 1 - start0 + FLANK)
        mask[lo:hi] = True
        got = {}
        for arm, seq in (("ref", ref_seq), ("alt", alt_seq)):
            o = aq.call_with_quota_retry(
                lambda s=seq: model.predict_sequence(sequence=s, requested_outputs=outs,
                                                    ontology_terms=LIVER_TERMS, interval=interval),
                label=f"{tag}:{arm}")
            got[arm] = summarise(o, mask)
        row = {"tag": tag, "chrom": chrom, "pos_hg38": pos1, "ref": ref, "alt": alt,
               "variant_class": indel_class(ref, alt), "length_delta": length_delta(ref, alt)}
        for ch in ("rna", "atac", "dnase", "h3k27ac", "splice"):
            a, b = got["alt"][ch], got["ref"][ch]
            row[f"{ch}_ref"] = b
            row[f"{ch}_alt"] = a
            row[f"{ch}_log2"] = math.log2((a + EPS) / (b + EPS)) if a == a and b == b else math.nan
        return row

    dbsnp = {}
    dp = sorted((la.out_root().parent).glob("p0-dbsnp-orientation-*/tables/dbsnp_indel_orientation.tsv"))
    dp = [q for q in dp if not (q.parent.parent / "SUPERSEDED.txt").exists()]
    if dp:
        dbsnp = dbsnp_orientations(la.read_tsv(dp[-1]))
        la.log(f"step 63: dbSNP orients {len(dbsnp)} deferred indels ({dp[-1]})")
    else:
        la.log("step 63: no dbSNP orientation table found; reference rule only")

    rows = []
    for sid in targets:
        c = xw.get(sid)
        parsed = parse_source_id(sid)
        chrom = (c or {}).get("hg38_chrom") or ""
        pos38 = (c or {}).get("hg38_position_1based") or ""
        if not chrom or not pos38:
            rows.append({"tag": sid, "chrom": parsed["chrom"], "pos_hg38": 0, "ref": parsed["ref"],
                         "alt": parsed["alt"], "variant_class": indel_class(parsed["ref"], parsed["alt"]),
                         "length_delta": length_delta(parsed["ref"], parsed["alt"]),
                         "state": "no_hg38_coordinate"})
            continue
        try:
            gdb = group_dbsnp_orientation(groups[sid]["members"], dbsnp)
            ref38, alt38, swap, amb, osrc = orientation_for(
                sid, chrom, int(pos38), parsed["ref"], parsed["alt"],
                lambda c, a, b: fasta.fetch(c, a, b), {sid: gdb} if gdb else {})
            r = score(chrom, int(pos38), ref38, alt38, sid)
            if r:
                r["allele_swap"] = swap
                r["orientation_ambiguous"] = amb
                r["orientation_source"] = osrc
                r["dbsnp_rule"] = gdb["dbsnp_rule"] if osrc == "dbsnp" else ""
        except la.ContractError as exc:
            rows.append({"tag": sid, "chrom": chrom, "pos_hg38": int(pos38), "ref": parsed["ref"],
                         "alt": parsed["alt"], "variant_class": indel_class(parsed["ref"], parsed["alt"]),
                         "length_delta": length_delta(parsed["ref"], parsed["alt"]),
                         "state": f"refused: {exc}"})
            continue
        if not r:
            continue
        gsigs = groups[sid]["signals"]
        r.update(state="scored", posterior_mass=groups[sid]["mass"], n_signals=len(gsigs),
                 source_variant_id=sid, arm="indel", member_source_ids=";".join(groups[sid]["members"]),
                 signals=";".join(sorted(gsigs)[:5]),
                 genes=";".join(sorted({signals[s]["gene"] for s in gsigs if s in signals})[:5]),
                 traits=";".join(sorted({signals[s]["trait"] for s in gsigs if s in signals})[:5]))
        rows.append(r)

        # matched SNV controls at the SAME signals, scored the same way: the magnitude is unreadable alone.
        pool = sorted({v for s in gsigs for _, v in snv_pool.get(s, [])})
        rng = np.random.default_rng(int(__import__("hashlib").sha256(sid.encode()).hexdigest()[:8], 16))
        for v in (list(rng.permutation(pool))[:N_SNV_CONTROLS] if pool else []):
            cc = xw.get(str(v))
            if not cc or not cc.get("hg38_chrom"):
                continue
            # A mapped SNV already carries its validated hg38 alleles in the crosswalk; re-parsing the
            # hg19 source id would reintroduce the same allele-order assumption that refused the 12.
            try:
                cr = score(cc["hg38_chrom"], int(cc["hg38_position_1based"]), cc["hg38_ref"], cc["hg38_alt"],
                           f"{sid}|ctrl|{v}")
                if cr:
                    cr["allele_swap"] = cc.get("allele_swap", "")
            except la.ContractError:
                continue
            if cr:
                cr.update(state="scored", arm="snv_control", control_for=sid, source_variant_id=str(v))
                rows.append(cr)

    cols = sorted({k for r in rows for k in r})
    la.write_tsv_once(TABLES / "indel_rescue_effects.tsv", rows, cols)

    def arr(arm, ch):
        return np.asarray([r[f"{ch}_log2"] for r in rows
                           if r.get("arm") == arm and isinstance(r.get(f"{ch}_log2"), float)
                           and r[f"{ch}_log2"] == r[f"{ch}_log2"]], dtype=float)

    summary = {
        "what_this_is": ("model-API predictions for indels the Atlas will not serve. These are predicted TRACK "
                         "ratios, not calibrated quantiles, and are never pooled with Atlas quantiles."),
        "atlas_cannot_reach_indels_by_either_route": ("point query returns UNIMPLEMENTED; interval saturation "
                                                      "returned 0 multi-base alleles in 104,352 variants over 40 "
                                                      "archives (checked 2026-09-13)"),
        "n_deferred_indel_ids_total": len(indel_mass),
        "n_distinct_deferred_indels": len(groups),
        "n_targets_above_mass_threshold": len(targets),
        "target_unit": ("the distinct variant: source ids naming one indel at one hg38 site with the same "
                        "unordered allele pair are merged, their mass summed and their signals pooled"),
        "min_posterior_mass": MIN_MASS,
        "n_scored": sum(r.get("state") == "scored" and r.get("arm") == "indel" for r in rows),
        "n_snv_controls": sum(r.get("arm") == "snv_control" for r in rows),
        "not_scored": [{"variant": r["tag"], "state": r["state"]} for r in rows
                       if r.get("state") not in (None, "scored")],
        "median_abs_log2": {ch: {arm: (float(np.median(np.abs(arr(arm, ch)))) if arr(arm, ch).size else None)
                                 for arm in ("indel", "snv_control")}
                            for ch in ("rna", "atac", "dnase", "h3k27ac", "splice")},
        "claim_boundary": ("a predicted effect in a reference context for one static sequence; it does not make the "
                           "variant causal, does not assign a gene, and cannot be compared with an Atlas quantile"),
    }
    json.dump(summary, (TABLES / "indel_rescue_summary.json").open("w"), indent=1, default=float)
    la.log(f"step 63: scored {summary['n_scored']} indels and {summary['n_snv_controls']} SNV controls")


if __name__ == "__main__":
    main()
