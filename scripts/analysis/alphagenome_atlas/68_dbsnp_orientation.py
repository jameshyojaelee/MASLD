#!/usr/bin/env python3
"""Step 68: resolve deferred-indel allele orientation against dbSNP.

The problem this closes. A VCF indel's short allele is a prefix of its long one, so when the long allele is
at the reference BOTH readings are representable and the reference cannot say which variant the GWAS tested:
100,468 of this Resource's 114,389 deferred indels (87.8%) are in that state. Step 67 fell back to the
source's allele order, which is a weak prior, and step 64 had to carry an `orientation_resolved` stratum of
only 16 paired targets as the honest subset.

dbSNP stores every variant in its canonical REF/ALT form, so a dbSNP record whose allele PAIR matches ours
names which of our two alleles is the reference. Matching is on CANONICAL form, not on strings: dbSNP collapses a repeat into one multi-allelic record whose
REF spans the whole run, so `C/CTT` sits inside rs3038218 as REF=CTTTTTTT with twenty ALTs. Step 69
normalises both sides (verified against `bcftools norm`). Without it 77% of the indels this step called
"absent from dbSNP" were simply written differently.

The primary arm is dbSNP build 157 on GRCh38.p14, normalised against this repository's single reference
FASTA and in the assembly the variants are scored in. Build 156 on GRCh37 (already on disk) confirms by
exact match only -- no hg19 FASTA is introduced, so it cannot left-shift through a repeat. A variant the two
assemblies disagree about is reported, never silently resolved.

dbSNP is an arbiter, not a repair: where it disagrees with the current reference-validated call, the
disagreement is written out rather than applied.

Outputs (tables/): dbsnp_indel_orientation.tsv, dbsnp_indel_orientation_summary.json
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import pathlib
import re
import sys
from collections import defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la


def _norm_mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("norm", pathlib.Path(__file__).resolve().parent / "69_normalise_indels.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


NORM = _norm_mod()

ROOT = la.out_root()
TABLES = ROOT / "tables"
TRACK0 = la.track0_root() / "tables"
SLIM_HG19 = pathlib.Path(os.environ.get(
    "AGA_DBSNP_HG19_SLIM",
    "/gpfs/commons/groups/sanjana_lab/mdrabkin/dbsnp_data/GCF_000001405.25_slim.tsv"))
VCF_HG38 = pathlib.Path(os.environ.get(
    "AGA_DBSNP_HG38_VCF",
    "/gpfs/commons/home/jameslee/reference_genome/dbsnp/GCF_000001405.40.gz"))

ACC_RE = re.compile(r"^NC_0000(\d{2})\.\d+$")
OUT_COLS = ["source_variant_id", "hg19_chrom", "hg19_pos", "source_allele1", "source_allele2",
            "hg38_chrom", "hg38_position_1based", "current_variant_uid", "current_orientation_ambiguous",
            "dbsnp_hg19_state", "dbsnp_hg19_ref", "dbsnp_hg19_alt", "dbsnp_hg19_rsid",
            "dbsnp_hg38_state", "dbsnp_hg38_ref", "dbsnp_hg38_alt", "dbsnp_hg38_rsid",
            "assemblies_agree", "verdict_vs_current", "resolved_variant_uid",
            "frequency_tie_break", "resolved_alt_frequency"]


FREQ_FLOOR = 0.01


def parse_freq(freq_field: str, n_alts: int) -> list:
    """Per-ALT allele frequency from dbSNP's FREQ, as `STUDY:refAF,alt1AF,...|STUDY2:...`.

    Returns one value per ALT (index 0 = first ALT), taking the maximum across studies and treating '.' as
    unreported. A study that lists fewer values than the record has ALTs is used for the alleles it does
    cover rather than being dropped or shifted.
    """
    out = [None] * int(n_alts)
    for study in str(freq_field or "").split("|"):
        if ":" not in study:
            continue
        values = study.split(":", 1)[1].split(",")
        for i, raw in enumerate(values[1:]):          # values[0] is the reference allele
            if i >= len(out) or raw in (".", ""):
                continue
            try:
                v = float(raw)
            except ValueError:
                continue
            out[i] = v if out[i] is None else max(out[i], v)
    return out


def tie_break_by_frequency(freq_by_orientation: dict, floor: float = FREQ_FLOOR) -> dict:
    """Pick the orientation whose alternate allele is common, where dbSNP holds both reciprocal events.

    In a homopolymer both the longer and shorter allele are real variants, so the allele list cannot say
    which one a GWAS tested. Frequency can: a GWAS indel is common by construction, so if exactly one of the
    two candidate alternate alleles clears the floor, that is the variant the study carried. A missing
    frequency is treated as not common, never as evidence for the allele.
    """
    common = [k for k, v in freq_by_orientation.items() if v is not None and float(v) >= floor]
    if len(common) == 1:
        ref, alt = common[0]
        return {"state": "resolved_by_frequency", "ref": ref, "alt": alt,
                "alt_frequency": float(freq_by_orientation[common[0]])}
    if len(common) > 1:
        return {"state": "both_common", "ref": None, "alt": None}
    return {"state": "neither_common", "ref": None, "alt": None}


def chrom_of(accession: str):
    """RefSeq chromosome accession -> chr name, or None for anything not a primary chromosome.

    NC_000001.10 is GRCh37 chr1 and NC_000001.11 is GRCh38 chr1: the number before the dot is the
    chromosome, the version is the assembly. Unplaced scaffolds are not chromosomes and get no name.
    """
    m = ACC_RE.match(str(accession).strip())
    if m is None:
        return None
    n = int(m.group(1))
    if 1 <= n <= 22:
        return f"chr{n}"
    if n == 23:
        return "chrX"
    if n == 24:
        return "chrY"
    return None


def orientation_from_records(records, a1: str, a2: str) -> dict:
    """Which of our two alleles dbSNP calls the reference, given the records at that position.

    `records` is an iterable of (REF, [ALT, ...]). A record matches when its REF is one of our alleles and
    the other is among its ALTs. Both orientations present means dbSNP holds two distinct variants there and
    the source id cannot be attributed to one of them; that is reported, not guessed.
    """
    x, y = str(a1).upper(), str(a2).upper()
    hits = set()
    for ref, alts in records:
        r = str(ref).upper()
        alt_set = {str(a).upper() for a in alts}
        if r == x and y in alt_set:
            hits.add((x, y))
        if r == y and x in alt_set:
            hits.add((y, x))
    if len(hits) == 1:
        ref, alt = hits.pop()
        return {"state": "resolved", "ref": ref, "alt": alt}
    if len(hits) > 1:
        return {"state": "both_in_dbsnp", "ref": None, "alt": None}
    return {"state": "absent", "ref": None, "alt": None}


def compare_to_current(current_uid: str, ref: str, alt: str) -> str:
    """Does dbSNP's orientation agree with the reference-validated call already written?

    Only the alleles are compared: the caller supplies a uid from the same crosswalk row, so the site is
    identical by construction and this function is never in a position to check it.
    """
    uid = str(current_uid or "").strip()
    if not uid:
        return "no_current_call"
    parts = uid.split(":")
    if len(parts) != 4:
        raise la.ContractError(f"variant uid is not chrom:pos:ref:alt: {uid!r}")
    cur_ref, cur_alt = parts[2].upper(), parts[3].upper()
    if {cur_ref, cur_alt} != {str(ref).upper(), str(alt).upper()}:
        return "different_alleles"
    return "agrees" if cur_ref == str(ref).upper() else "flips"


def reconcile(o19: dict, o38: dict, hg38_checked: bool) -> tuple:
    """Combine the two assembly arms into one orientation, or refuse.

    hg38 is the arbiter: build 157, normalised against this repository's single reference FASTA, with the
    frequency tie-break. hg19 confirms by exact match only. A non-answer from either arm is NOT a
    disagreement -- only two resolved answers that differ are, and those are refused.
    """
    r19, r38 = o19["state"] == "resolved", o38["state"] == "resolved"
    if not hg38_checked:
        return "not_checked", (o19 if r19 else None)
    if r19 and r38:
        if (o19["ref"], o19["alt"]) == (o38["ref"], o38["alt"]):
            return "agree", o38
        return "DISAGREE", None
    if r38:
        return "one_resolved_only", o38
    if r19:
        return "one_resolved_only", o19
    return "neither_resolved", None


def unresolved_label(o19: dict, o38: dict) -> str:
    """Name the state that actually blocked resolution, preferring the arbiter's."""
    if o19["state"] == "resolved" and o38["state"] == "resolved":
        return "refused_assemblies_disagree"
    state = o38["state"] if o38["state"] != "resolved" else o19["state"]
    return "unresolved_" + state


def load_targets() -> dict:
    """Deferred indels, keyed by source id, with their hg19 site, alleles and current hg38 call."""
    out = {}
    with la.open_text(TRACK0 / "variant_crosswalk.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["exclusion_reason"] not in ("indel_deferred", "indel_deferred_ref_mismatch"):
                continue
            m = la.SOURCE_ID_RE.match(r["source_variant_id"])
            if m is None:
                raise la.ContractError(f"unparsable source id {r['source_variant_id']!r}")
            a1, a2 = r["source_alleles"].split("/")
            out[r["source_variant_id"]] = {
                "hg19_chrom": f"chr{m['chrom']}", "hg19_pos": int(m["pos"]),
                "source_allele1": a1.upper(), "source_allele2": a2.upper(),
                "hg38_chrom": r["hg38_chrom"], "hg38_position_1based": r["hg38_position_1based"]}
    return out


def scan_slim_hg19(sites: set) -> dict:
    """One streaming pass over the GRCh37 CHROM/POS/ID/REF/ALT table, keeping only our positions."""
    got = defaultdict(list)
    with SLIM_HG19.open() as h:
        rd = csv.reader(h, delimiter="\t")
        next(rd, None)
        for row in rd:
            if len(row) < 5:
                continue
            c = chrom_of(row[0])
            if c is None:
                continue
            key = (c, int(row[1]))
            if key in sites:
                got[key].append((int(row[1]), row[3], row[4].split(","), row[2]))
    return got


def query_vcf_hg38(sites_by_chrom: dict) -> dict:
    """Position queries against the GRCh38 dbSNP VCF via its tabix index."""
    import pysam

    got = defaultdict(list)
    vcf = pysam.VariantFile(str(VCF_HG38))
    contigs = {chrom_of(c): c for c in vcf.header.contigs if chrom_of(c)}
    for chrom, positions in sorted(sites_by_chrom.items()):
        acc = contigs.get(chrom)
        if acc is None:
            continue
        for pos in sorted(positions):
            try:
                # A repeat record's REF can START before our position and still describe our event, so the
                # window is widened and normalisation decides, not the start coordinate.
                for rec in vcf.fetch(acc, max(0, pos - 60), pos + 60):
                    alts = list(rec.alts or [])
                    raw = rec.info.get("FREQ")
                    field = ",".join(raw) if isinstance(raw, (tuple, list)) else (raw or "")
                    got[(chrom, pos)].append((rec.pos, rec.ref, alts, rec.id or "",
                                              parse_freq(field, len(alts))))
            except ValueError:
                continue
    return got


# One reference genome for every sequence operation in this repository. No hg19 FASTA is introduced: the
# hg19 arm therefore matches representations exactly (no left-shift through a repeat) and is a secondary
# confirmation, while hg38 -- build 157, verified base-for-base against this FASTA -- is the primary arm.
HG38_FASTA = ("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
              "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")


def main() -> None:
    import pysam

    fa38 = pysam.FastaFile(HG38_FASTA)
    def fetch38(c, s, e):
        return fa38.fetch(c, s, e)
    def fetch19(c, s, e):
        return ""          # no hg19 FASTA: trimming only, never a left-shift

    targets = load_targets()
    la.log(f"step 68: {len(targets)} deferred indels to orient")
    sites19 = {(t["hg19_chrom"], t["hg19_pos"]) for t in targets.values()}
    if not SLIM_HG19.exists():
        raise la.ContractError(f"GRCh37 dbSNP table not found: {SLIM_HG19}")
    hg19 = scan_slim_hg19(sites19)
    la.log(f"step 68: hg19 dbSNP records at {len(hg19)} of {len(sites19)} positions")

    hg38 = {}
    if VCF_HG38.exists() and (VCF_HG38.parent / (VCF_HG38.name + ".tbi")).exists():
        by_chrom = defaultdict(set)
        for t in targets.values():
            if t["hg38_chrom"] and t["hg38_position_1based"]:
                by_chrom[t["hg38_chrom"]].add(int(t["hg38_position_1based"]))
        hg38 = query_vcf_hg38(by_chrom)
        la.log(f"step 68: hg38 dbSNP records at {len(hg38)} positions")
    else:
        la.log(f"step 68: hg38 dbSNP VCF absent ({VCF_HG38}); hg19 pass only, no cross-assembly confirmation")

    repair = {}
    rp = sorted((ROOT.parent).glob("p0-indel-uid-repair-*/tables/deferred_indel_uid_repair.tsv"))
    rp = [q for q in rp if not (q.parent.parent / "SUPERSEDED.txt").exists()]
    if rp:
        for r in la.read_tsv(rp[-1]):
            repair[r["source_variant_id"]] = r

    rows = []
    for sid, t in sorted(targets.items()):
        cur = repair.get(sid, {})
        cur_uid = cur.get("recovered_variant_uid", "")
        raw19 = hg19.get((t["hg19_chrom"], t["hg19_pos"]), [])
        ids19 = [i for _, _, _, i in raw19]
        o19 = NORM.orientation_from_normalised([(p, r, a) for p, r, a, _ in raw19], t["hg19_chrom"],
                                               t["hg19_pos"], t["source_allele1"], t["source_allele2"],
                                               fetch19)
        key38 = (t["hg38_chrom"], int(t["hg38_position_1based"])) if t["hg38_position_1based"] else None
        raw38 = hg38.get(key38, []) if key38 else []
        ids38 = [i for _, _, _, i, _ in raw38]
        o38 = (NORM.orientation_from_normalised([(p, r, a) for p, r, a, _, _ in raw38], t["hg38_chrom"],
                                                int(t["hg38_position_1based"]), t["source_allele1"],
                                                t["source_allele2"], fetch38)
               if key38 else {"state": "absent", "ref": None, "alt": None})
        # Where dbSNP holds both reciprocal events, allele frequency says which one a GWAS could carry.
        freq_state, freq_alt = "", ""
        if o38["state"] == "both_in_dbsnp":
            fb = NORM.frequency_by_orientation([(p, r, a) for p, r, a, _, _ in raw38], t["hg38_chrom"],
                                               int(t["hg38_position_1based"]), t["source_allele1"],
                                               t["source_allele2"], fetch38,
                                               [f for _, _, _, _, f in raw38])
            tb = tie_break_by_frequency(fb)
            freq_state = tb["state"]
            if tb["state"] == "resolved_by_frequency":
                o38 = {"state": "resolved", "ref": tb["ref"], "alt": tb["alt"]}
                freq_alt = repr(tb["alt_frequency"])

        agree, chosen = reconcile(o19, o38, bool(hg38))
        verdict = (compare_to_current(cur_uid, chosen["ref"], chosen["alt"]) if chosen
                   else unresolved_label(o19, o38))
        resolved_uid = ""
        if chosen and t["hg38_chrom"] and t["hg38_position_1based"]:
            resolved_uid = la.variant_uid(t["hg38_chrom"], int(t["hg38_position_1based"]),
                                          chosen["ref"], chosen["alt"])
        rows.append({"source_variant_id": sid, **{k: t[k] for k in
                     ("hg19_chrom", "hg19_pos", "source_allele1", "source_allele2",
                      "hg38_chrom", "hg38_position_1based")},
                     "current_variant_uid": cur_uid,
                     "current_orientation_ambiguous": cur.get("orientation_ambiguous", ""),
                     "dbsnp_hg19_state": o19["state"], "dbsnp_hg19_ref": o19["ref"] or "",
                     "dbsnp_hg19_alt": o19["alt"] or "", "dbsnp_hg19_rsid": ";".join(sorted(set(ids19))[:3]),
                     "dbsnp_hg38_state": o38["state"], "dbsnp_hg38_ref": o38["ref"] or "",
                     "dbsnp_hg38_alt": o38["alt"] or "", "dbsnp_hg38_rsid": ";".join(sorted(set(ids38))[:3]),
                     "assemblies_agree": agree, "verdict_vs_current": verdict,
                     "resolved_variant_uid": resolved_uid,
                     "frequency_tie_break": freq_state, "resolved_alt_frequency": freq_alt})

    la.write_tsv_once(TABLES / "dbsnp_indel_orientation.tsv", rows, OUT_COLS)
    amb = [r for r in rows if r["current_orientation_ambiguous"] == "True"]
    res = [r for r in rows if r["resolved_variant_uid"]]
    summary = {
        "n_deferred_indels": len(rows),
        "hg19_source": str(SLIM_HG19), "hg38_source": str(VCF_HG38) if hg38 else None,
        "dbsnp_hg19_states": {s: sum(r["dbsnp_hg19_state"] == s for r in rows)
                              for s in sorted({r["dbsnp_hg19_state"] for r in rows})},
        "dbsnp_hg38_states": {s: sum(r["dbsnp_hg38_state"] == s for r in rows)
                              for s in sorted({r["dbsnp_hg38_state"] for r in rows})},
        "assemblies_agree": {s: sum(r["assemblies_agree"] == s for r in rows)
                             for s in sorted({r["assemblies_agree"] for r in rows})},
        "verdict_vs_current": {s: sum(r["verdict_vs_current"] == s for r in rows)
                               for s in sorted({r["verdict_vs_current"] for r in rows})},
        "n_resolved_by_dbsnp": len(res),
        "frequency_tie_break": {s: sum(r["frequency_tie_break"] == s for r in rows)
                                for s in sorted({r["frequency_tie_break"] for r in rows}) if s},
        "frequency_floor": FREQ_FLOOR,
        "why_frequency_can_arbitrate": ("in a homopolymer both the longer and shorter allele are real "
                                        "variants, so the allele list cannot say which one a GWAS tested; "
                                        "a GWAS indel is common by construction, so a single allele clearing "
                                        "the floor identifies it. gnomAD is not needed: dbSNP ships these "
                                        "frequencies in its FREQ field."),
        "previously_ambiguous": {
            "n": len(amb),
            "now_resolved": sum(1 for r in amb if r["resolved_variant_uid"]),
            "dbsnp_flips_the_current_call": sum(1 for r in amb if r["verdict_vs_current"] == "flips")},
        "dbsnp_is_an_arbiter_not_a_repair": ("rows where dbSNP disagrees with the reference-validated call "
                                             "are written out with verdict_vs_current == 'flips'; nothing is "
                                             "overwritten here"),
    }
    json.dump(summary, (TABLES / "dbsnp_indel_orientation_summary.json").open("w"), indent=1, default=float)
    la.log(f"step 68: resolved {len(res)} of {len(rows)}; of the {len(amb)} previously ambiguous, "
           f"{summary['previously_ambiguous']['now_resolved']} resolved, "
           f"{summary['previously_ambiguous']['dbsnp_flips_the_current_call']} flip the current call")


if __name__ == "__main__":
    main()
