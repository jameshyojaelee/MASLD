#!/usr/bin/env python3
"""Step 67: recover the hg38 identity of the deferred indels in an already-written Track 0 run.

The crosswalk kept indels out of the Atlas query set, which is correct — the point-query API refuses them —
but it also returned before writing their hg38 alleles, so `variant_uid` is empty for all 114,389 of them.
A posterior keyed on that column collapses every deferred indel of a signal into one entry whose weight is
their sum. `lib_atlas.crosswalk_variant` now fills the uid, so any future run is correct; this step repairs
the run already on disk, which is write-once and will not be regenerated.

Every one of those 114,389 rows already has a 1:1 liftover and an hg38 coordinate. Only the alleles were
missing, and they are recovered fail-closed against the reference FASTA: a row whose reference allele is not
actually at that position keeps no uid and is recorded as `ref_mismatch`.

Recovering an identity does NOT make a variant queryable. `atlas_queryable` stays False on every repaired
row; what changes is that the variant can be counted, joined to the model-API rescue, and distinguished
from the other deferred indels at its signal.

Outputs (tables/): deferred_indel_uid_repair.tsv, deferred_indel_uid_repair_summary.json
"""

from __future__ import annotations

import csv
import json
import os
import pathlib

import pysam

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
TRACK0 = la.track0_root() / "tables"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
OUT_COLS = ["source_variant_id", "variant_uid", "recovered_variant_uid", "identity_after_repair", "repair_state",
            "hg38_chrom", "hg38_position_1based", "hg38_ref", "hg38_alt", "allele_swap", "variant_class",
            "mapping_status", "exclusion_reason", "atlas_queryable", "orientation_ambiguous", "n_signals"]


def source_alleles(field: str) -> tuple:
    parts = str(field).split("/")
    if len(parts) != 2 or not all(parts):
        raise la.ContractError(f"source_alleles is not a pair: {field!r}")
    return parts[0], parts[1]


def _variant_class(ref: str, alt: str) -> str:
    r, a = str(ref).upper(), str(alt).upper()
    if len(r) == len(a):
        return "snv" if len(r) == 1 else "mnv"
    return "insertion" if len(a) > len(r) else "deletion"


def repair_row(row: dict, fetch) -> dict:
    """One crosswalk row, with its identity recovered where the reference supports it."""
    out = {"source_variant_id": row["source_variant_id"], "variant_uid": row.get("variant_uid", "") or "",
           "recovered_variant_uid": "", "repair_state": "", "hg38_chrom": row.get("hg38_chrom", "") or "",
           "hg38_position_1based": row.get("hg38_position_1based", "") or "", "hg38_ref": "", "hg38_alt": "",
           "allele_swap": "", "variant_class": "", "orientation_ambiguous": "",
           "mapping_status": row.get("mapping_status", ""),
           "exclusion_reason": row.get("exclusion_reason", "") or "",
           "atlas_queryable": str(row.get("mapping_status") == "mapped"),
           "n_signals": row.get("n_signals", "")}
    if out["variant_uid"]:
        a1, a2 = source_alleles(row["source_alleles"])
        out.update(repair_state="already_identified", recovered_variant_uid=out["variant_uid"],
                   variant_class=_variant_class(*out["variant_uid"].split(":")[2:4]))
        return out
    if not (out["hg38_chrom"] and out["hg38_position_1based"]):
        out["repair_state"] = "no_hg38_coordinate"
        return out
    if str(row.get("n_liftover_mappings", "")) != "1":
        out["repair_state"] = "not_one_to_one"
        return out
    a1, a2 = source_alleles(row["source_alleles"])
    r = la.resolve_indel_alleles(out["hg38_chrom"], int(out["hg38_position_1based"]), a1, a2, fetch)
    if r["variant_uid"] is None:
        out["repair_state"] = "ref_mismatch"
        return out
    out.update(recovered_variant_uid=r["variant_uid"], hg38_ref=r["hg38_ref"], hg38_alt=r["hg38_alt"],
               allele_swap=r["allele_swap"], variant_class=_variant_class(r["hg38_ref"], r["hg38_alt"]),
               orientation_ambiguous=str(bool(r["orientation_ambiguous"])), repair_state="recovered")
    return out


def assert_identities_agree(rows) -> int:
    """Rows sharing an identity must describe the same hg38 site; returns how many identities are shared.

    The crosswalk is keyed on the SOURCE string, and universe A spells a chromosome `10` where universe B
    spells it `chr10`, so one variant reached through both universes legitimately has two rows. A signal
    belongs to exactly one universe, so no posterior can contain both spellings. Two rows that disagree on
    the site would be a different thing entirely and abort the run.
    """
    from collections import defaultdict

    groups = defaultdict(list)
    for r in rows:
        groups[r["identity_after_repair"]].append(r)
    shared = 0
    for key, grp in groups.items():
        if len(grp) == 1:
            continue
        shared += 1
        sites = {(g["hg38_chrom"], str(g["hg38_position_1based"]), g["hg38_ref"], g["hg38_alt"]) for g in grp}
        if len(sites) > 1:
            raise la.ContractError(f"identity {key} names {len(sites)} different hg38 sites: {sorted(sites)}")
    return shared


def signal_identities(weight_rows, repair: dict) -> list:
    """(signal, identity) pairs a posterior would use once the repair is applied."""
    out = []
    for r in weight_rows:
        uid = (r.get("variant_uid") or "").strip() or repair.get(r["source_variant_id"], "")
        out.append((r["signal_uid"], la.posterior_key(uid, r["source_variant_id"])))
    return out


def assert_no_merge(pairs) -> None:
    """The repair must not reintroduce the bug it fixes.

    A variant written with its alleles in both orders (`T:TG` and `TG:T` at one hg19 position) resolves to a
    single reference-oriented uid, which is the point. If both spellings ever appeared in ONE signal's
    posterior, keying on that uid would sum two entries into one, exactly as the empty uid did.
    """
    la.assert_unique_identities(pairs)


def identity_after_repair(row: dict) -> str:
    """The key a posterior should use once the repair is applied."""
    return la.posterior_key(row.get("variant_uid") or row.get("recovered_variant_uid") or "",
                            row.get("source_variant_id", ""))


def main() -> None:
    fasta = pysam.FastaFile(FASTA)
    only = os.environ.get("AGA_REPAIR_ONLY_UNIDENTIFIED", "1") == "1"
    rows, states = [], {}
    with la.open_text(TRACK0 / "variant_crosswalk.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if only and (r.get("variant_uid") or "").strip():
                states["already_identified"] = states.get("already_identified", 0) + 1
                continue
            out = repair_row(r, fasta.fetch)
            out["identity_after_repair"] = identity_after_repair(out)
            states[out["repair_state"]] = states.get(out["repair_state"], 0) + 1
            rows.append(out)

    n_shared = assert_identities_agree(rows)
    # The guard above is a consistency check on the writer and cannot fail for a recovered row, whose
    # identity IS its site. The guard that matters reads the posterior itself.
    repair = {r["source_variant_id"]: r["recovered_variant_uid"] for r in rows if r["repair_state"] == "recovered"}
    with la.open_text(TRACK0 / "signal_variant_weights.tsv.gz") as h:
        assert_no_merge(signal_identities(csv.DictReader(h, delimiter="\t"), repair))
    if any(r["atlas_queryable"] != "False" for r in rows):
        raise la.ContractError("a repaired row was marked queryable; recovering an identity does not make "
                               "the point-query API accept an indel")
    la.write_tsv_once(TABLES / "deferred_indel_uid_repair.tsv", rows, OUT_COLS)

    recovered = [r for r in rows if r["repair_state"] == "recovered"]
    summary = {
        "repairs_run_against": str(TRACK0 / "variant_crosswalk.tsv.gz"),
        "source_sha256": la.sha256_file(TRACK0 / "variant_crosswalk.tsv.gz"),
        "n_rows_examined": len(rows), "repair_states": states,
        "n_identities_shared_by_two_source_spellings": n_shared,
        "why_an_identity_can_be_shared": ("the crosswalk is keyed on the source string and universe A spells a "
                                          "chromosome `10` where universe B spells it `chr10`, so a variant "
                                          "reached through both universes has two rows; they are required to "
                                          "agree on the hg38 site and a signal belongs to only one universe"),
        "n_recovered": len(recovered),
        "recovered_by_class": {c: sum(r["variant_class"] == c for r in recovered)
                               for c in sorted({r["variant_class"] for r in recovered})},
        "n_orientation_ambiguous": sum(r["orientation_ambiguous"] == "True" for r in recovered),
        "n_orientation_resolved_by_the_reference": sum(r["orientation_ambiguous"] == "False" for r in recovered),
        "what_ambiguous_means": ("a VCF indel's short allele is a prefix of its long one, so when the long "
                                 "allele is at the reference both readings are representable. The reported "
                                 "reference allele then rests on the SOURCE's allele order, not on the "
                                 "reference, and the two readings give different alternate sequences."),
        "what_this_does_not_change": ("the Atlas query set, any served prediction, or any adopted number. "
                                      "Steps 08/09/34/62 already filtered to mapped variants or fell back to "
                                      "the source id, so no written table keyed on the empty uid."),
        "what_it_does_change": ("every deferred indel now has a distinct hg38 identity, so a posterior can be "
                               "keyed on it without merging, and the model-API rescue joins on a coordinate "
                               "rather than on a normalised hg19 string"),
        "atlas_queryable_after_repair": sorted({r["atlas_queryable"] for r in rows}),
        "no_signal_merges_two_source_variants_onto_one_recovered_uid": True,
    }
    json.dump(summary, (TABLES / "deferred_indel_uid_repair_summary.json").open("w"), indent=1, default=float)
    la.log(f"step 67: {len(rows)} rows examined, {len(recovered)} identities recovered; states {states}")


if __name__ == "__main__":
    main()
