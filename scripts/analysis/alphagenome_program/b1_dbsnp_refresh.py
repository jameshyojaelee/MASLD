#!/usr/bin/env python3
"""B1 dbSNP refresh: recompute the research layer's indel-orientation columns against the CURRENT,
non-superseded dbSNP orientation deposit, and state what moves.

Reads only. Nothing under GWAS/finemapping/results/alphagenome_atlas/ is written or edited.

The four dbSNP orientation deposits, oldest first:
  p0-dbsnp-orientation-20260914T121856Z  SUPERSEDED (string REF/ALT match; 24,902 indels wrongly absent)
  p0-dbsnp-orientation-20260914T174542Z  SUPERSEDED (no frequency tie-break; 68,296 left both_in_dbsnp)
  p0-dbsnp-orientation-20260914T192341Z  SUPERSEDED (hg19 non-answer treated as DISAGREE; vetoed 66,859)
  p0-dbsnp-orientation-20260914T210732Z  CURRENT, no SUPERSEDED.txt

The v2 research layer b1-response-layer-v2-20260914T204700Z deposited variant_class /
orientation_source / dbsnp_deposit built on the UNION of the two oldest (earlier then later wins).
This script:
  1. reproduces those deposited columns under the union with an independent loader (guard);
  2. recomputes them under each deposit alone and under the current deposit;
  3. counts the rows that change, for all 10,511 layer rows and for the 548 indel rows;
  4. does the same for the P6f rescue targets and the P6d dossier uids;
  5. pins the HSD17B13 rs72613567 anchor across all four deposits.

Usage: b1_dbsnp_refresh.py <output directory>
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import pathlib
import re
import subprocess
import sys

csv.field_size_limit(10_000_000)

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS = PROJECT / "GWAS/finemapping/results/alphagenome_atlas"
PROGRAM = PROJECT / "GWAS/finemapping/results/alphagenome_program"

V2_LAYER = PROGRAM / "b1-response-layer-v2-20260914T204700Z"
V2_VARIANTS = V2_LAYER / "response_layer_variants.tsv.gz"
V2_PROVENANCE = V2_LAYER / "indel_orientation_provenance.tsv"

DEPOSITS = {
    "earlier_20260914T121856Z": ATLAS / "p0-dbsnp-orientation-20260914T121856Z",
    "later_20260914T174542Z": ATLAS / "p0-dbsnp-orientation-20260914T174542Z",
    "rerun_20260914T192341Z": ATLAS / "p0-dbsnp-orientation-20260914T192341Z",
    "current_20260914T210732Z": ATLAS / "p0-dbsnp-orientation-20260914T210732Z",
}
CURRENT = "current_20260914T210732Z"
UNION_PARTS = ("earlier_20260914T121856Z", "later_20260914T174542Z")  # later wins where it speaks

REPAIR_TSV = ATLAS / "p0-indel-uid-repair-20260913T193750Z/tables/deferred_indel_uid_repair.tsv"
RESCUE_152 = ATLAS / "p6f-indel-rescue-20260914T135052Z/tables/indel_rescue_effects.tsv"
RESCUE_NEW = ATLAS / "p6f-indel-rescue-20260914T231617Z/tables/indel_rescue_effects.tsv"
DOSSIER_USED_BY_LAYER = ATLAS / "p6d-dossier-20260914T173952Z/tables/clinical_variant_dossier.tsv"
DOSSIER_NEWEST = ATLAS / "p6d-dossier-20260915T000218Z/tables/clinical_variant_dossier.tsv"

# The deposits spell most source ids without "chr" (73 of 114,389 carry it), P6f and P6d mostly with it.
# Every lookup here is raw-first then normalised; the normalisation is checked for collisions that would
# merge two different calls before it is used.
HSD_SOURCE = "hg19:4:88231392:T:TA"
HSD_EXPECTED_UID = "chr4:87310240:T:TA"

SOURCE_ID_RE = re.compile(
    r"^(?P<asm>[A-Za-z0-9]+):(?:chr)?(?P<chrom>[0-9XYMxym]+|MT):(?P<pos>\d+):(?P<ref>[A-Za-z]+):(?P<alt>[A-Za-z]+)$")


class ContractError(RuntimeError):
    pass


# ---------------------------------------------------------------- identity helpers
# Copied verbatim from b1_response_layer_v2.py so the refresh uses the SAME loader semantics as the
# deposited columns. The reproduction guard below fails loudly if they have drifted.

def normalise_source_id(source_id: str) -> str:
    m = SOURCE_ID_RE.match(str(source_id).strip())
    if m is None:
        raise ContractError(f"source id is not assembly:chrom:pos:ref:alt: {source_id!r}")
    return f"{m['asm']}:chr{m['chrom']}:{int(m['pos'])}:{m['ref'].upper()}:{m['alt'].upper()}"


def variant_class(ref: str, alt: str) -> str:
    r, a = str(ref).upper(), str(alt).upper()
    if len(r) == len(a):
        return "snv" if len(r) == 1 else "mnv"
    return "insertion" if len(a) > len(r) else "deletion"


def orientation_source(dbsnp_uid: str, recovered_uid: str) -> str:
    if str(dbsnp_uid or "").strip():
        return "dbsnp"
    if str(recovered_uid or "").strip():
        return "reference_or_source_order"
    return ""


def alleles_of(variant_uid, source_variant_id, recovered_uid="", dbsnp_uid=""):
    for uid in (str(variant_uid or "").strip(), str(dbsnp_uid or "").strip(),
                str(recovered_uid or "").strip()):
        if uid:
            parts = uid.split(":")
            if len(parts) != 4:
                raise ContractError(f"variant uid is not chrom:pos:ref:alt: {uid!r}")
            return parts[2], parts[3]
    parts = normalise_source_id(source_variant_id).split(":")
    return parts[3], parts[4]


def resolved_identity(variant_uid, source_variant_id, recovered_uid="", dbsnp_uid=""):
    """The hg38 identity the layer would print for this row, in the layer's own precedence order."""
    for uid in (str(variant_uid or "").strip(), str(dbsnp_uid or "").strip(),
                str(recovered_uid or "").strip()):
        if uid:
            return uid
    return ""


# ---------------------------------------------------------------- io

def read_tsv(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv(path, rows, columns=None):
    columns = columns or list(rows[0].keys())
    with open(path, "w") as h:
        w = csv.DictWriter(h, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in columns})
    return len(rows)


def sha256(path, cap=400_000_000):
    p = pathlib.Path(path)
    if not p.exists() or p.stat().st_size > cap:
        return ""
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------- deposits

class Lookup:
    """Raw-first, normalised-fallback map from a source id to a deposit row.

    The source id is spelled two ways across these deposits (`hg19:4:...` and `hg19:chr4:...`). A
    normalised-only map would merge the two spellings; a raw-only map would drop every P6f / P6d join.
    Collisions where the two spellings give DIFFERENT calls are counted and refused.
    """

    def __init__(self, rows, value_of):
        self.raw, byn, self.collisions, self.conflicts = {}, {}, 0, []
        for r in rows:
            v = value_of(r)
            self.raw[r["source_variant_id"]] = v
            k = normalise_source_id(r["source_variant_id"])
            if k in byn:
                self.collisions += 1
                if byn[k] != v:
                    self.conflicts.append(k)
            else:
                byn[k] = v
        self.norm = byn

    def get(self, src, default=""):
        if src in self.raw:
            return self.raw[src]
        try:
            return self.norm.get(normalise_source_id(src), default)
        except ContractError:
            return default

    def __contains__(self, src):
        return bool(self.get(src, ""))

    def __len__(self):
        return len(self.norm)


class UnionLookup:
    """The layer's deposited rule: the later deposit wins where it speaks, the earlier fills in."""

    def __init__(self, earlier, later, earlier_name, later_name):
        self.earlier, self.later = earlier, later
        self.earlier_name, self.later_name = earlier_name, later_name

    def get(self, src, default=""):
        v = self.later.get(src, "")
        if v:
            return v
        return self.earlier.get(src, "") or default

    def which(self, src):
        if self.later.get(src, ""):
            return self.later_name
        return self.earlier_name if self.earlier.get(src, "") else ""

    def __contains__(self, src):
        return bool(self.get(src, ""))

    def __len__(self):
        return len(set(self.earlier.norm) | set(self.later.norm))


def load_deposits():
    out = {}
    for name, d in DEPOSITS.items():
        tsv = d / "tables/dbsnp_indel_orientation.tsv"
        if not tsv.exists() or tsv.stat().st_size == 0:
            raise ContractError(f"{tsv} absent or empty")
        rows = read_tsv(tsv)
        resolved = [r for r in rows if r["resolved_variant_uid"]]
        uid = Lookup(resolved, lambda r: r["resolved_variant_uid"])
        detail = Lookup(rows, lambda r: r)
        if uid.conflicts:
            raise ContractError(f"{name}: {len(uid.conflicts)} source ids resolve differently under the "
                                f"two chr spellings, e.g. {uid.conflicts[0]}")
        summ = json.load((d / "tables/dbsnp_indel_orientation_summary.json").open())
        out[name] = {
            "dir": d, "superseded": (d / "SUPERSEDED.txt").exists(),
            "superseded_text": ((d / "SUPERSEDED.txt").read_text().strip()
                                if (d / "SUPERSEDED.txt").exists() else ""),
            "n_rows": len(rows), "n_resolved": len(resolved), "uid": uid, "detail": detail,
            "rows": rows, "spelling_collisions": uid.collisions, "summary": summ,
        }
        log(f"  {name}: {len(rows)} rows, {len(resolved)} resolved, "
            f"superseded={out[name]['superseded']}, chr-spelling collisions={uid.collisions} "
            f"(0 conflicting)")
    return out


SEED = 20260915
N_BOOT = 2000


def flip_rate_with_blocks(dep_rows, label):
    """The flip rate against the source allele order, resampled on the 1-Mb block, not the row.

    The earlier deposit's 38,438 flips against 38,342 agrees was read as an exact coin flip. Source ids
    inside one locus are LD-correlated and are not independent draws, so the unit is the 1-Mb hg38 block.
    Reported as the row-level rate, the block-mean rate, and a seeded block bootstrap interval on the
    block-mean.
    """
    import random
    blocks = {}
    n_flip = n_agree = 0
    for r in dep_rows:
        v = r["verdict_vs_current"]
        if v not in ("flips", "agrees"):
            continue
        chrom = r.get("hg38_chrom") or ""
        pos = r.get("hg38_position_1based") or ""
        if not chrom or not pos:
            continue
        b = f"{chrom}~{int(pos) // 1_000_000}"
        f, a = blocks.get(b, (0, 0))
        blocks[b] = (f + (v == "flips"), a + (v == "agrees"))
        n_flip += (v == "flips")
        n_agree += (v == "agrees")
    keys = sorted(blocks)
    per_block = [blocks[k][0] / (blocks[k][0] + blocks[k][1]) for k in keys]
    block_mean = sum(per_block) / len(per_block) if per_block else float("nan")
    rng = random.Random(SEED)
    draws = []
    for _ in range(N_BOOT):
        s = [per_block[rng.randrange(len(per_block))] for _ in range(len(per_block))]
        draws.append(sum(s) / len(s))
    draws.sort()
    lo, hi = draws[int(0.025 * N_BOOT)], draws[int(0.975 * N_BOOT) - 1]
    return {
        "deposit": label,
        "n_rows_with_a_verdict": n_flip + n_agree,
        "n_flips": n_flip, "n_agrees": n_agree,
        "row_level_flip_rate": round(n_flip / (n_flip + n_agree), 6) if (n_flip + n_agree) else "",
        "n_1mb_blocks": len(keys),
        "block_mean_flip_rate": round(block_mean, 6),
        "block_bootstrap_95ci_low": round(lo, 6), "block_bootstrap_95ci_high": round(hi, 6),
        "n_boot": N_BOOT, "seed": SEED,
        "half_inside_the_interval": str(lo <= 0.5 <= hi),
    }


def deposit_comparison(dep):
    """Deposit-level resolution, side by side. n_deferred_indels is the same denominator in all four."""
    rows = []
    for name in ["earlier_20260914T121856Z", "later_20260914T174542Z",
                 "rerun_20260914T192341Z", CURRENT]:
        d = dep[name]
        s = d["summary"]
        n_def = s["n_deferred_indels"]
        rows.append({
            "deposit": name,
            "state": ("SUPERSEDED" if d["superseded"] else "CURRENT"),
            "superseded_reason": d["superseded_text"].replace("\n", " "),
            "n_deferred_indels": n_def,
            "n_resolved_by_dbsnp": s["n_resolved_by_dbsnp"],
            "resolution_rate": round(s["n_resolved_by_dbsnp"] / n_def, 6),
            "n_resolved_rows_in_table": d["n_resolved"],
            "verdict_agrees": s["verdict_vs_current"]["agrees"],
            "verdict_flips": s["verdict_vs_current"]["flips"],
            "flip_rate_among_resolved": round(
                s["verdict_vs_current"]["flips"]
                / (s["verdict_vs_current"]["flips"] + s["verdict_vs_current"]["agrees"]), 6),
            "unresolved_absent": s["verdict_vs_current"]["unresolved_absent"],
            "unresolved_both_in_dbsnp": s["verdict_vs_current"]["unresolved_both_in_dbsnp"],
            "previously_ambiguous_now_resolved": s["previously_ambiguous"]["now_resolved"],
            "hg38_state_resolved": s["dbsnp_hg38_states"]["resolved"],
            "hg38_state_both_in_dbsnp": s["dbsnp_hg38_states"]["both_in_dbsnp"],
            "hg38_state_absent": s["dbsnp_hg38_states"]["absent"],
            "frequency_resolved": s.get("frequency_tie_break", {}).get("resolved_by_frequency", ""),
        })
    return rows


# ---------------------------------------------------------------- research layer recompute

ARMS = ["earlier_20260914T121856Z", "later_20260914T174542Z", "union_deposited",
        "rerun_20260914T192341Z", CURRENT]


def recompute_layer(dep, repair):
    union = UnionLookup(dep[UNION_PARTS[0]]["uid"], dep[UNION_PARTS[1]]["uid"],
                        "dbsnp_20260914T121856Z", "dbsnp_20260914T174542Z")
    maps = {a: (union if a == "union_deposited" else dep[a]["uid"]) for a in ARMS}

    layer = read_tsv(V2_VARIANTS)
    out = []
    for r in layer:
        uid, src = r["variant_uid"], r["source_variant_id"]
        rec = repair.get(src, "")
        row = {"signal_uid": r["signal_uid"], "variant_key": r["variant_key"],
               "source_variant_id": src, "variant_uid": uid, "hg38_uid_recovered": rec,
               "deposited_variant_class": r["variant_class"],
               "deposited_orientation_source": r["orientation_source"],
               "deposited_dbsnp_deposit": r["dbsnp_deposit"],
               "deposited_hg38_uid_dbsnp": r["hg38_uid_dbsnp"],
               "deposited_hg38_ref": r["hg38_ref"], "deposited_hg38_alt": r["hg38_alt"]}
        for a in ARMS:
            dbu = maps[a].get(src, "")
            ref, alt = alleles_of(uid, src, rec, dbu)
            row[f"uid_dbsnp__{a}"] = dbu
            row[f"resolved_identity__{a}"] = resolved_identity(uid, src, rec, dbu)
            row[f"variant_class__{a}"] = variant_class(ref, alt)
            row[f"orientation_source__{a}"] = orientation_source(dbu, rec)
            row[f"hg38_ref__{a}"] = ref
            row[f"hg38_alt__{a}"] = alt
        row["dbsnp_deposit__union_deposited"] = union.which(src)
        out.append(row)
    return out, union


def reproduction_guard(rows):
    """The deposited columns must come back out of this loader under the union arm, on every row."""
    bad = {"variant_class": [], "orientation_source": [], "dbsnp_deposit": [],
           "hg38_uid_dbsnp": [], "hg38_ref": [], "hg38_alt": []}
    for r in rows:
        if r["variant_class__union_deposited"] != r["deposited_variant_class"]:
            bad["variant_class"].append(r["variant_key"])
        if r["orientation_source__union_deposited"] != r["deposited_orientation_source"]:
            bad["orientation_source"].append(r["variant_key"])
        if r["dbsnp_deposit__union_deposited"] != r["deposited_dbsnp_deposit"]:
            bad["dbsnp_deposit"].append(r["variant_key"])
        if r["uid_dbsnp__union_deposited"] != r["deposited_hg38_uid_dbsnp"]:
            bad["hg38_uid_dbsnp"].append(r["variant_key"])
        if r["hg38_ref__union_deposited"] != r["deposited_hg38_ref"]:
            bad["hg38_ref"].append(r["variant_key"])
        if r["hg38_alt__union_deposited"] != r["deposited_hg38_alt"]:
            bad["hg38_alt"].append(r["variant_key"])
    return {"n_rows": len(rows),
            "n_mismatched": {k: len(v) for k, v in bad.items()},
            "examples": {k: v[:3] for k, v in bad.items() if v},
            "verdict": ("PASS" if not any(bad.values()) else "FAIL")}


def same_site_reciprocal(a: str, b: str) -> bool:
    """Both uids name the same chrom and position with the allele pair swapped: the reciprocal variant.

    This is the change that reverses the sign of every signed readout, because the effect was scored on
    ALT-over-REF of the other allele order.
    """
    if not a or not b or a == b:
        return False
    pa, pb = a.split(":"), b.split(":")
    if len(pa) != 4 or len(pb) != 4:
        return False
    return pa[0] == pb[0] and pa[1] == pb[1] and pa[2] == pb[3] and pa[3] == pb[2]


def change_counts(rows, ref_arm, label):
    """How many rows move between ref_arm and the current deposit, on each of the three properties."""
    n_cls = n_ident = n_orient = n_uid = 0
    n_recip = n_other_move = 0
    cls_transitions, ident_examples = {}, []
    for r in rows:
        if r[f"variant_class__{ref_arm}"] != r[f"variant_class__{CURRENT}"]:
            n_cls += 1
            k = f'{r[f"variant_class__{ref_arm}"]}->{r[f"variant_class__{CURRENT}"]}'
            cls_transitions[k] = cls_transitions.get(k, 0) + 1
        a, b = r[f"resolved_identity__{ref_arm}"], r[f"resolved_identity__{CURRENT}"]
        if a != b:
            n_ident += 1
            if same_site_reciprocal(a, b):
                n_recip += 1
            else:
                n_other_move += 1
            if len(ident_examples) < 5:
                ident_examples.append({"variant_key": r["variant_key"], ref_arm: a, CURRENT: b,
                                       "same_site_reciprocal": same_site_reciprocal(a, b)})
        if r[f"orientation_source__{ref_arm}"] != r[f"orientation_source__{CURRENT}"]:
            n_orient += 1
        if r[f"uid_dbsnp__{ref_arm}"] != r[f"uid_dbsnp__{CURRENT}"]:
            n_uid += 1
    return {"comparison": label, "n_rows": len(rows),
            "n_variant_class_changes": n_cls,
            "variant_class_transitions": cls_transitions,
            "n_resolved_identity_changes": n_ident,
            "n_identity_changes_that_are_the_same_site_reciprocal": n_recip,
            "n_identity_changes_that_move_position_or_alleles_otherwise": n_other_move,
            "n_orientation_source_changes": n_orient,
            "n_dbsnp_uid_changes": n_uid,
            "resolved_identity_examples": ident_examples}


# ---------------------------------------------------------------- P6f

def p6f_targets(path, dep, label):
    rows = read_tsv(path)
    ind = [r for r in rows if r.get("arm") == "indel"]
    cur = dep[CURRENT]["uid"]
    curdet = dep[CURRENT]["detail"]
    out, n_join, n_change, n_class_change, n_newly_oriented = [], 0, 0, 0, 0
    distinct = set()
    for r in ind:
        src = r["source_variant_id"]
        used_uid = f'{r["chrom"]}:{r["pos_hg38"]}:{r["ref"]}:{r["alt"]}'
        distinct.add(used_uid)
        # a merged target carries every source id that named it; any of them can carry the dbSNP call
        members = [m for m in (r.get("member_source_ids") or "").split(";") if m] or [src]
        cur_uid = ""
        for m in [src] + members:
            cur_uid = cur.get(m, "")
            if cur_uid:
                break
        joined = bool(cur_uid)
        n_join += joined
        changed = joined and cur_uid != used_uid
        n_change += changed
        cls_used = variant_class(r["ref"], r["alt"])
        cls_cur = variant_class(*cur_uid.split(":")[2:4]) if cur_uid else cls_used
        n_class_change += (joined and cls_cur != cls_used)
        if joined and r.get("orientation_source") != "dbsnp":
            n_newly_oriented += 1
        recip = same_site_reciprocal(used_uid, cur_uid)
        det = {}
        for m in [src] + members:
            det = curdet.get(m, {}) or {}
            if det:
                break
        out.append({
            "deposit": label, "source_variant_id": src,
            "member_source_ids": r.get("member_source_ids", ""),
            "uid_used_by_p6f": used_uid,
            "orientation_source_in_p6f": r.get("orientation_source", ""),
            "orientation_ambiguous_in_p6f": r.get("orientation_ambiguous", ""),
            "variant_class_in_p6f": cls_used,
            "uid_under_current_deposit": cur_uid,
            "variant_class_under_current_deposit": cls_cur,
            "current_deposit_joins": str(joined),
            "orientation_changes_under_current": str(bool(changed)),
            "same_site_reciprocal_under_current": str(recip),
            "variant_class_changes_under_current": str(bool(joined and cls_cur != cls_used)),
            "current_dbsnp_hg38_state": det.get("dbsnp_hg38_state", ""),
            "current_frequency_tie_break": det.get("frequency_tie_break", ""),
            "current_resolved_alt_frequency": det.get("resolved_alt_frequency", ""),
            "current_verdict_vs_reference_call": det.get("verdict_vs_current", ""),
            "splice_log2": r.get("splice_log2", ""), "rna_log2": r.get("rna_log2", ""),
            "posterior_mass": r.get("posterior_mass", ""), "genes": r.get("genes", ""),
        })
    summary = {"deposit": label, "n_indel_rows": len(ind),
               "n_distinct_uids_as_scored": len(distinct),
               "n_distinct_source_ids_normalised": len({normalise_source_id(r["source_variant_id"])
                                                        for r in ind}),
               "n_joining_current_deposit": n_join,
               "n_orientation_changes_under_current": n_change,
               "n_of_those_that_are_the_same_site_reciprocal":
                   sum(1 for r in out if r["same_site_reciprocal_under_current"] == "True"),
               "n_variant_class_changes_under_current": n_class_change,
               "n_not_previously_dbsnp_oriented_that_current_resolves": n_newly_oriented,
               "n_distinct_uids_changing": len({r["uid_used_by_p6f"] for r in out
                                                if r["orientation_changes_under_current"] == "True"})}
    return out, summary


# ---------------------------------------------------------------- P6d

def p6d_uids(path, dep, repair, label):
    rows = read_tsv(path)
    cur = dep[CURRENT]["uid"]
    union = UnionLookup(dep[UNION_PARTS[0]]["uid"], dep[UNION_PARTS[1]]["uid"],
                        "dbsnp_20260914T121856Z", "dbsnp_20260914T174542Z")
    out = []
    n_dbsnp_used = n_change_uid = n_change_class = n_newly = 0
    for r in rows:
        src, uid = r["source_variant_id"], r["variant_uid"]
        rec = repair.get(src, r.get("hg38_uid_recovered", ""))
        used_db = r.get("hg38_uid_dbsnp", "")
        cur_db = cur.get(src, "")
        used_ident = resolved_identity(uid, src, rec, used_db)
        cur_ident = resolved_identity(uid, src, rec, cur_db)
        used_cls = variant_class(*alleles_of(uid, src, rec, used_db))
        cur_cls = variant_class(*alleles_of(uid, src, rec, cur_db))
        if used_db:
            n_dbsnp_used += 1
            if cur_db and cur_db != used_db:
                n_change_uid += 1
        if not used_db and cur_db:
            n_newly += 1
        if cur_ident != used_ident:
            pass
        if cur_cls != used_cls:
            n_change_class += 1
        out.append({
            "deposit": label, "signal_uid": r["signal_uid"], "source_variant_id": src,
            "variant_uid": uid, "hg38_uid_recovered": rec,
            "hg38_uid_dbsnp_in_dossier": used_db,
            "uid_dbsnp_under_current_deposit": cur_db,
            "uid_dbsnp_under_union": union.get(src, ""),
            "resolved_identity_in_dossier": used_ident,
            "resolved_identity_under_current": cur_ident,
            "variant_class_in_dossier": r.get("variant_class", ""),
            "variant_class_recomputed_in_dossier": used_cls,
            "variant_class_under_current": cur_cls,
            "orientation_source_in_dossier": r.get("orientation_source", ""),
            "orientation_source_under_current": orientation_source(cur_db, rec),
            "identity_changes_under_current": str(cur_ident != used_ident),
            "same_site_reciprocal_under_current": str(same_site_reciprocal(used_ident, cur_ident)),
            "variant_class_changes_under_current": str(cur_cls != used_cls),
        })
    n_ident_change = sum(1 for r in out if r["identity_changes_under_current"] == "True")
    summary = {"deposit": label, "n_rows": len(rows),
               "n_distinct_uids_in_dossier": len({r["resolved_identity_in_dossier"] for r in out
                                                  if r["resolved_identity_in_dossier"]}),
               "n_rows_with_dbsnp_uid_in_dossier": n_dbsnp_used,
               "n_distinct_dbsnp_uids_in_dossier": len({r["hg38_uid_dbsnp_in_dossier"] for r in out
                                                        if r["hg38_uid_dbsnp_in_dossier"]}),
               "n_of_those_that_change_under_current": n_change_uid,
               "n_distinct_dbsnp_uids_that_change": len({r["hg38_uid_dbsnp_in_dossier"] for r in out
                                                         if r["hg38_uid_dbsnp_in_dossier"]
                                                         and r["uid_dbsnp_under_current_deposit"]
                                                         and r["uid_dbsnp_under_current_deposit"]
                                                         != r["hg38_uid_dbsnp_in_dossier"]}),
               "n_rows_current_resolves_that_dossier_did_not": n_newly,
               "n_resolved_identity_changes": n_ident_change,
               "n_of_those_that_are_the_same_site_reciprocal":
                   sum(1 for r in out if r["same_site_reciprocal_under_current"] == "True"),
               "n_variant_class_changes": n_change_class}
    return out, summary


# ---------------------------------------------------------------- anchor

def anchor_check(dep):
    rows = []
    for name in ["earlier_20260914T121856Z", "later_20260914T174542Z",
                 "rerun_20260914T192341Z", CURRENT]:
        d = dep[name]["detail"].get(HSD_SOURCE, {})
        rows.append({
            "deposit": name,
            "state": ("SUPERSEDED" if dep[name]["superseded"] else "CURRENT"),
            "source_variant_id": HSD_SOURCE,
            "present_in_deposit": str(bool(d)),
            "resolved_variant_uid": d.get("resolved_variant_uid", ""),
            "matches_expected_anchor": str(d.get("resolved_variant_uid", "") == HSD_EXPECTED_UID),
            "dbsnp_hg38_state": d.get("dbsnp_hg38_state", ""),
            "dbsnp_hg38_rsid": d.get("dbsnp_hg38_rsid", ""),
            "assemblies_agree": d.get("assemblies_agree", ""),
            "verdict_vs_current": d.get("verdict_vs_current", ""),
            "frequency_tie_break": d.get("frequency_tie_break", ""),
            "resolved_alt_frequency": d.get("resolved_alt_frequency", ""),
        })
    return rows


# ---------------------------------------------------------------- main

def main():
    out = pathlib.Path(sys.argv[1]).resolve()
    tables = out / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    log("== loading dbSNP deposits")
    dep = load_deposits()

    log("== deposit-level comparison")
    cmp_rows = deposit_comparison(dep)
    write_tsv(tables / "dbsnp_deposit_resolution_comparison.tsv", cmp_rows)
    for r in cmp_rows:
        log(f"  {r['deposit']:<28} {r['state']:<11} resolved {r['n_resolved_by_dbsnp']:>7} "
            f"({r['resolution_rate']:.4f})  flips {r['verdict_flips']:>6} "
            f"(rate {r['flip_rate_among_resolved']:.4f})")

    log("== flip rate against the source allele order, resampled on the 1-Mb block")
    flip_rows = [flip_rate_with_blocks(dep[n]["rows"], n)
                 for n in ["earlier_20260914T121856Z", "later_20260914T174542Z",
                           "rerun_20260914T192341Z", CURRENT]]
    write_tsv(tables / "flip_rate_by_1mb_block.tsv", flip_rows)
    for f in flip_rows:
        log(f"  {f['deposit']:<28} rows {f['n_rows_with_a_verdict']:>6} flip {f['row_level_flip_rate']:.4f} "
            f"| {f['n_1mb_blocks']} blocks, block-mean {f['block_mean_flip_rate']:.4f} "
            f"[{f['block_bootstrap_95ci_low']:.4f}, {f['block_bootstrap_95ci_high']:.4f}] "
            f"0.5 inside={f['half_inside_the_interval']}")

    repair = Lookup([r for r in read_tsv(REPAIR_TSV) if r["repair_state"] == "recovered"],
                    lambda r: r["recovered_variant_uid"])
    if repair.conflicts:
        raise ContractError(f"the repair table gives {len(repair.conflicts)} source ids two identities "
                            f"under the two chr spellings, e.g. {repair.conflicts[0]}")
    log(f"== deferred-indel repair: {len(repair)} recovered identities "
        f"({repair.collisions} chr-spelling collisions, 0 conflicting)")

    log("== recomputing the research layer under every deposit")
    rows, union = recompute_layer(dep, repair)
    log(f"  layer rows: {len(rows)}; union resolves {len(union)} source ids; "
        f"current resolves {len(dep[CURRENT]['uid'])}")

    guard = reproduction_guard(rows)
    log(f"== reproduction guard (union arm vs deposited columns): {guard['verdict']} "
        f"{guard['n_mismatched']}")
    if guard["verdict"] != "PASS":
        raise ContractError(f"the refresh loader does not reproduce the deposited columns: {guard}")

    indel = [r for r in rows if r["deposited_variant_class"] in ("insertion", "deletion")]
    log(f"  indel rows in the layer (deposited class): {len(indel)}")
    prov = read_tsv(V2_PROVENANCE)
    prov_keys = {(r["signal_uid"], r["variant_key"]) for r in prov}
    mine_keys = {(r["signal_uid"], r["variant_key"]) for r in indel}
    log(f"  deposited provenance rows: {len(prov)}; key sets identical: {prov_keys == mine_keys}")

    # a row could in principle enter or leave the indel set under the current deposit
    indel_current = [r for r in rows if r[f"variant_class__{CURRENT}"] in ("insertion", "deletion")]
    log(f"  indel rows under the current deposit: {len(indel_current)}")

    # A row that already carries an Atlas-served variant_uid takes its alleles from that uid, so the dbSNP
    # call cannot move its class or identity however much the arbiter changes. The rows that CAN move are
    # the deferred ones, whose uid is empty.
    deferred = [r for r in indel if not r["variant_uid"]]
    mapped = [r for r in indel if r["variant_uid"]]
    log(f"  of the {len(indel)} indel rows: {len(deferred)} deferred (empty variant_uid, the dbSNP call "
        f"can move them) and {len(mapped)} already carrying an hg38 uid (it cannot)")

    changes = []
    for arm, label in [("earlier_20260914T121856Z", "current vs earlier 121856Z alone"),
                       ("later_20260914T174542Z", "current vs later 174542Z alone"),
                       ("union_deposited", "current vs the deposited union (earlier + later)"),
                       ("rerun_20260914T192341Z", "current vs rerun 192341Z alone")]:
        changes.append({"scope": "all_layer_rows", **change_counts(rows, arm, label)})
        changes.append({"scope": "indel_rows_548", **change_counts(indel, arm, label)})
        changes.append({"scope": "indel_rows_deferred_only", **change_counts(deferred, arm, label)})
    for c in changes:
        log(f"  [{c['scope']}] {c['comparison']}: class {c['n_variant_class_changes']}, "
            f"identity {c['n_resolved_identity_changes']}, "
            f"orientation_source {c['n_orientation_source_changes']}, "
            f"dbsnp uid {c['n_dbsnp_uid_changes']}")

    # refreshed column set: the indel rows, all arms side by side
    cols = ["signal_uid", "variant_key", "source_variant_id", "variant_uid", "hg38_uid_recovered",
            "deposited_variant_class", "deposited_orientation_source", "deposited_dbsnp_deposit",
            "deposited_hg38_uid_dbsnp", "deposited_hg38_ref", "deposited_hg38_alt"]
    for a in ARMS:
        cols += [f"uid_dbsnp__{a}", f"resolved_identity__{a}", f"variant_class__{a}",
                 f"orientation_source__{a}", f"hg38_ref__{a}", f"hg38_alt__{a}"]
    cols += ["dbsnp_deposit__union_deposited"]
    write_tsv(tables / "research_layer_indel_orientation_refreshed.tsv", indel, cols)
    write_tsv(tables / "research_layer_orientation_refreshed_all_rows.tsv", rows, cols)

    log("== P6f rescue targets")
    p6f_rows, p6f_sum = p6f_targets(RESCUE_152, dep, "p6f-indel-rescue-20260914T135052Z")
    p6f_rows_new, p6f_sum_new = p6f_targets(RESCUE_NEW, dep, "p6f-indel-rescue-20260914T231617Z")
    write_tsv(tables / "p6f_orientation_under_current_deposit.tsv", p6f_rows + p6f_rows_new)
    log(f"  {json.dumps(p6f_sum)}")
    log(f"  {json.dumps(p6f_sum_new)}")

    log("== P6d dossier uids")
    p6d_rows, p6d_sum = p6d_uids(DOSSIER_USED_BY_LAYER, dep, repair, "p6d-dossier-20260914T173952Z")
    p6d_rows_new, p6d_sum_new = p6d_uids(DOSSIER_NEWEST, dep, repair, "p6d-dossier-20260915T000218Z")
    write_tsv(tables / "p6d_orientation_under_current_deposit.tsv", p6d_rows + p6d_rows_new)
    log(f"  {json.dumps(p6d_sum)}")
    log(f"  {json.dumps(p6d_sum_new)}")

    log("== HSD17B13 rs72613567 anchor")
    anc = anchor_check(dep)
    write_tsv(tables / "hsd17b13_anchor_across_deposits.tsv", anc)
    for a in anc:
        log(f"  {a['deposit']:<28} {a['resolved_variant_uid'] or '(absent)':<24} "
            f"matches={a['matches_expected_anchor']} state={a['dbsnp_hg38_state']}")

    summary = {
        "package": "b1-dbsnp-refresh",
        "current_deposit": DEPOSITS[CURRENT].name,
        "current_deposit_is_superseded": dep[CURRENT]["superseded"],
        "superseded_deposits": {k: dep[k]["superseded_text"] for k in DEPOSITS if dep[k]["superseded"]},
        "deposit_comparison": cmp_rows,
        "flip_rate_by_1mb_block": flip_rows,
        "research_layer": {
            "layer": V2_LAYER.name,
            "n_rows": len(rows),
            "n_indel_rows_deposited_class": len(indel),
            "n_indel_rows_under_current_deposit": len(indel_current),
            "n_indel_rows_deferred_empty_variant_uid": len(deferred),
            "n_indel_rows_already_carrying_an_hg38_uid": len(mapped),
            "class_counts_per_arm": {
                a: {c: sum(1 for r in indel if r[f"variant_class__{a}"] == c)
                    for c in sorted({r[f"variant_class__{a}"] for r in indel})} for a in ARMS},
            "n_indel_rows_resolved_per_arm": {
                a: sum(1 for r in indel if r[f"uid_dbsnp__{a}"]) for a in ARMS},
            "provenance_key_sets_identical": prov_keys == mine_keys,
            "reproduction_guard": guard,
            "changes": changes,
            "n_source_ids_the_union_resolves": len(union),
            "n_source_ids_the_current_deposit_resolves": len(dep[CURRENT]["uid"]),
            "chr_spelling_collisions_per_deposit": {k: dep[k]["spelling_collisions"] for k in DEPOSITS},
        },
        "p6f": [p6f_sum, p6f_sum_new],
        "p6d": [p6d_sum, p6d_sum_new],
        "hsd17b13_anchor": anc,
        "what_this_does_not_do": (
            "no deposit under alphagenome_atlas/ is edited; the v2 layer's own columns are not rewritten; "
            "no model API call is made and no effect size is rescored. The refreshed columns say which "
            "allele each row's uid names under the current arbiter, not whether an effect is real."),
    }
    json.dump(summary, (tables / "b1_dbsnp_refresh_summary.json").open("w"), indent=1)

    # manifest
    man = []
    for p in [V2_VARIANTS, V2_PROVENANCE, REPAIR_TSV, RESCUE_152, RESCUE_NEW,
              DOSSIER_USED_BY_LAYER, DOSSIER_NEWEST,
              pathlib.Path(__file__).resolve()] + \
             [d / "tables/dbsnp_indel_orientation.tsv" for d in DEPOSITS.values()] + \
             [d / "tables/dbsnp_indel_orientation_summary.json" for d in DEPOSITS.values()]:
        man.append({"path": str(p), "role": "input", "exists": str(p.exists()),
                    "bytes": (p.stat().st_size if p.exists() else 0), "sha256": sha256(p)})
    for p in sorted(tables.glob("*")):
        man.append({"path": str(p), "role": "output", "exists": "True",
                    "bytes": p.stat().st_size, "sha256": sha256(p)})
    write_tsv(out / "MANIFEST.tsv", man, ["path", "role", "exists", "bytes", "sha256"])

    (out / "pip_freeze.txt").write_text(
        subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True,
                       text=True).stdout)
    (out / "python_version.txt").write_text(sys.version + "\n" + sys.executable + "\n")
    log("== done")


if __name__ == "__main__":
    main()
