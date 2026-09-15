#!/usr/bin/env python3
"""Step 04: hg38 variant crosswalk, signal weights, analysis blocks, query set, pilot blocks."""

from __future__ import annotations

import csv
import json
from collections import defaultdict

import pysam

import lib_atlas as la

P = la.prespec()
OUT = la.out_root() / "tables"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"

XW_COLS = ["variant_uid", "source_assembly", "source_variant_id", "source_alleles", "hg38_chrom", "hg38_position_1based",
           "hg38_ref", "hg38_alt", "orientation", "allele_swap", "palindromic", "is_snv", "n_liftover_mappings",
           "mapping_status", "exclusion_reason", "n_signals", "posterior_key"]
W_COLS = ["signal_uid", "universe", "variant_uid", "posterior_key", "source_variant_id", "weight", "mapping_status",
          "exclusion_reason", "in_query_set"]


def main() -> None:
    fasta = pysam.FastaFile(FASTA)
    fetch = fasta.fetch
    signals = {s["signal_uid"]: s for s in la.read_tsv(OUT / "eligible_signals_provisional.tsv")}
    crosswalk: dict[str, dict] = {}          # keyed by source key -> crosswalk row
    weights: dict[str, dict[str, float]] = defaultdict(dict)   # signal -> source key -> weight

    def xw(source_key, assembly, chrom, pos, a1, a2, n_map):
        if source_key in crosswalk:
            return crosswalk[source_key]
        if chrom is None or pos in (None, "", "NA"):
            row = {"hg38_chrom": "", "hg38_position_1based": "", "source_alleles": f"{a1}/{a2}", "n_liftover_mappings": n_map,
                   "hg38_ref": None, "hg38_alt": None, "orientation": None, "allele_swap": None, "palindromic": False,
                   "is_snv": len(a1) == 1 and len(a2) == 1, "variant_uid": None, "mapping_status": "excluded",
                   "exclusion_reason": "liftover_failed" if n_map == 0 else "multimapped"}
        else:
            row = la.crosswalk_variant(chrom, int(pos), a1, a2, fetch, n_map)
        row.update(source_assembly=assembly, source_variant_id=source_key, n_signals=0)
        # Never blank. A variant the Atlas defers or fails to lift has no hg38 uid, and keying a posterior
        # on that column alone merges every one of them into a single summed entry per signal.
        row["posterior_key"] = la.posterior_key(row["variant_uid"] or "", source_key)
        crosswalk[source_key] = row
        return row

    # Universe A/C rows (already hg38 in the audit; re-validated against the FASTA here).
    with la.open_text(OUT / "variant_posteriors_A.tsv.gz") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            key = f"hg19:{r['hg19_variant_id']}:{r['allele1']}:{r['allele2']}"
            n_map = 1 if r["mapping_status"] == "mapped" else 0
            row = xw(key, "hg19", r["hg38_chrom"] or None, r["hg38_position_1based"], r["allele1"], r["allele2"], n_map)
            weights[r["signal_uid"]][key] = float(r["snp_pp_h4"])
    # Universe B rows.
    lift = {r["variant_key"]: r for r in la.read_tsv(OUT / "universeB_liftover_raw.tsv")}
    for r in la.read_tsv(OUT / "variant_posteriors_B_hg19.tsv"):
        vk = f"chr{r['chromosome']}:{r['position']}:{r['effect_allele']}:{r['other_allele']}"
        m = lift[vk]
        n_map = int(m["n_liftover_mappings"])
        chrom = m["hg38_chrom"] if n_map == 1 else None
        row = xw(f"hg19:{vk}", "hg19", chrom, m["hg38_position_1based"] if n_map == 1 else None,
                 r["effect_allele"], r["other_allele"], n_map)
        weights[r["signal_uid"]][f"hg19:{vk}"] = float(r["pip"])

    # Signal weights keyed by variant_uid; mass assertions; query floor.
    weight_rows, cov_rows, anchors = [], [], []
    query: set[str] = set()
    q = P["query_floor"]
    for sig, wmap in weights.items():
        universe = signals[sig]["universe"]
        mapped = {k: w for k, w in wmap.items() if crosswalk[k]["mapping_status"] == "mapped"}
        sel, excluded_floor = la.select_query_variants(mapped, q["weight_floor"], q["cumulative_mass"], q["per_signal_cap"])
        la.assert_signal_mass([(sig, crosswalk[k]["variant_uid"], w) for k, w in mapped.items()], tolerance=1e-6 if universe != "B_direct" else 0.1)
        for k, w in wmap.items():
            c = crosswalk[k]
            c["n_signals"] += 1
            weight_rows.append({"signal_uid": sig, "universe": universe, "variant_uid": c["variant_uid"] or "",
                                "posterior_key": c["posterior_key"], "source_variant_id": k, "weight": w,
                                "mapping_status": c["mapping_status"],
                                "exclusion_reason": c["exclusion_reason"] or "", "in_query_set": k in sel})
        for k in sel:
            query.add(crosswalk[k]["variant_uid"])
        total = sum(wmap.values())
        cov_rows.append({"signal_uid": sig, "universe": universe, "n_variants": len(wmap), "total_mass": total,
                         "mapped_mass": sum(mapped.values()), "unmapped_mass": total - sum(mapped.values()),
                         "n_queried": len(sel), "queried_mass": sum(mapped[k] for k in sel),
                         "excluded_mass_below_floor": excluded_floor})
        if mapped:
            top = max(mapped.items(), key=lambda kv: (kv[1], kv[0]))[0]
            anchors.append((sig, crosswalk[top]["hg38_chrom"], crosswalk[top]["hg38_position_1based"]))
    la.assert_unique_identities([(r["signal_uid"], r["posterior_key"]) for r in weight_rows])
    blocks = la.coarse_blocks(anchors)

    # Tables.
    xw_rows = sorted(crosswalk.values(), key=lambda c: (c["hg38_chrom"] or "zz", int(c["hg38_position_1based"] or 0), c["source_variant_id"]))
    la.write_tsv_once(OUT / "variant_crosswalk.tsv.gz", xw_rows, XW_COLS)
    la.write_tsv_once(OUT / "signal_variant_weights.tsv.gz", weight_rows, W_COLS)
    la.write_tsv_once(OUT / "signal_query_coverage.tsv", cov_rows, list(cov_rows[0].keys()))
    sig_rows = []
    for sig, s in signals.items():
        s = dict(s)
        s["analysis_block"] = blocks.get(sig, "")
        anchor = next(((c, p) for x, c, p in anchors if x == sig), ("", ""))
        s["anchor_hg38_chrom"], s["anchor_hg38_position"] = anchor
        sig_rows.append(s)
    la.write_tsv_once(OUT / "eligible_signals.tsv", sig_rows, list(sig_rows[0].keys()))
    block_rows = defaultdict(lambda: {"n_signals": 0, "universes": set(), "start": 10**12, "end": 0, "chrom": ""})
    for sig, chrom, pos in anchors:
        b = block_rows[blocks[sig]]
        b["n_signals"] += 1; b["universes"].add(signals[sig]["universe"]); b["chrom"] = chrom
        b["start"] = min(b["start"], int(pos)); b["end"] = max(b["end"], int(pos))
    ordered = sorted(block_rows.items(), key=lambda kv: (la.STANDARD_CHROMS and int(kv[1]["chrom"][3:]) if kv[1]["chrom"][3:].isdigit() else 99, kv[1]["start"]))
    la.write_tsv_once(OUT / "analysis_blocks.tsv",
                      [{"block_uid": k, "chrom": v["chrom"], "start": v["start"], "end": v["end"], "n_signals": v["n_signals"],
                        "universes": ";".join(sorted(v["universes"]))} for k, v in ordered],
                      ["block_uid", "chrom", "start", "end", "n_signals", "universes"])
    qrows = [{"variant_uid": u, "hg38_chrom": u.split(":")[0], "hg38_position_1based": int(u.split(":")[1]),
              "hg38_ref": u.split(":")[2], "hg38_alt": u.split(":")[3]} for u in sorted(query)]
    la.write_tsv_once(OUT / "query_set.tsv", qrows, ["variant_uid", "hg38_chrom", "hg38_position_1based", "hg38_ref", "hg38_alt"])
    direct_blocks = [k for k, v in ordered if v["universes"] & {"A_direct", "B_direct"}][:5]
    pilot_sigs = [s for s in sig_rows if s["analysis_block"] in direct_blocks]
    pilot_variants = sorted({w["variant_uid"] for w in weight_rows if w["signal_uid"] in {s["signal_uid"] for s in pilot_sigs} and w["in_query_set"]})
    with open(OUT / "pilot_blocks.json", "w") as handle:
        json.dump({"blocks": direct_blocks, "signals": [s["signal_uid"] for s in pilot_sigs], "variants": pilot_variants}, handle, indent=1)
    excl = defaultdict(int)
    for c in crosswalk.values():
        excl[c["exclusion_reason"] or "mapped"] += 1
    la.log(f"crosswalk: {len(crosswalk)} unique source variants; {dict(excl)}; query set {len(query)} unique hg38 SNVs; "
           f"blocks {len(block_rows)}; pilot blocks {direct_blocks} with {len(pilot_variants)} variants")


if __name__ == "__main__":
    main()
