#!/usr/bin/env python3
"""Step 42 (P4): full saturation retrieval over the frozen region universe.

Every possible substitution in every region of the universe fixed by step 41, retrieved with
`query_interval` and archived per region. Resumable: a region whose archive already carries a request.json
is skipped, so a restart continues rather than repeating.

Retrieval is NOT ontology-filtered server-side (a 50-curie `ontology_terms` list makes query_interval fail
with UNAVAILABLE; see the execution log). The filter is applied at archive time instead, which keeps only
the liver and proxy panel the analysis reads and takes a region from 46.7 MB to 11.0 MB. Track counts before
and after go into each request.json so the filter is auditable from the deposit.

Sized from the step-40 pilot: median 3.18 s and 11.0 MB per region, so about 91 hours and 1.14 TB over the
103,390-region universe.

usage: 42_saturation_retrieve.py [max_regions]
Outputs: raw/saturation/<universe>/<chrom_start-end>/, tables/saturation_retrieval_progress.json
"""

from __future__ import annotations

import json
import sys
import time

from alphagenome.data import genome

import atlas_archive as aa
import atlas_query as aq
import lib_atlas as la

ROOT = la.out_root()
RAW = ROOT / "raw"
TABLES = ROOT / "tables"
UNIVERSE = TABLES / "saturation_universe.tsv.gz"
SCORERS = ["ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE", "AVI_SCORE"]
ONTOLOGY = {
    "liver": ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182", "CL:0000632", "CL:0002571"],
    "HepG2": ["EFO:0001187"],
    "fibroblast_proxy": ["CL:0002547", "CL:0002548", "CL:0002549", "CL:0002550", "CL:0002551", "CL:0002552",
                          "CL:0002553", "CL:0002554", "CL:0002555", "CL:0002556", "CL:0002557", "CL:0002558"],
    "monocyte_macrophage_proxy": ["CL:0000235", "CL:0000576", "CL:0000862", "CL:0000863", "CL:0001054", "CL:0002057", "CL:0002397"],
    "endothelial_proxy": ["CL:0000115", "CL:0002138", "CL:0002144", "CL:0002188", "CL:0002543", "CL:0002544",
                           "CL:0002545", "CL:0002618", "CL:1000413", "CL:1000892", "CL:1001568", "CL:2000010"],
    "t_nk_proxy": ["CL:0000084", "CL:0000097", "CL:0000136", "CL:0000351", "CL:0000623", "CL:0000624",
                    "CL:0000625", "CL:0000792", "CL:0000798", "CL:0000823", "CL:0000895", "CL:0000897"],
}
TERMS = sorted({c for v in ONTOLOGY.values() for c in v})
REPORT_EVERY = 100


def load_universe() -> list:
    """Regions in the order step 41 fixed them, deduplicated on (universe, chrom, start, end)."""
    seen, out = set(), []
    for r in la.read_tsv(UNIVERSE):
        key = (r["universe"], r["chrom"], int(r["start0"]), int(r["end"]))
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def main() -> None:
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    regions = load_universe()
    if limit:
        regions = regions[:limit]
    client, sdk = aq.create_client()
    la.log(f"saturation: {len(regions)} regions in the frozen universe")
    done = skipped = 0
    t_start = time.time()
    for i, (uni, chrom, start0, end) in enumerate(regions):
        key = f"{chrom}:{start0}-{end}"
        d = RAW / "saturation" / uni / key.replace(":", "_")
        if (d / "request.json").exists():
            skipped += 1
            continue
        interval = genome.Interval(chromosome=chrom, start=start0, end=end)
        t0 = time.time()
        res = aq.call_with_quota_retry(
            lambda: client.query_interval(interval, requested_scorers=SCORERS, progress_bar=False, max_workers=4),
            label=key)
        elapsed = time.time() - t0
        res, track_counts = aq.filter_tracks_to_ontology(res, TERMS)
        aa.archive_scores(res, d, {"universe": uni, "interval": key, "requested_scorers": SCORERS,
                                   "ontology_terms_requested": None, "ontology_filter": "archive-time, client-side",
                                   "ontology_terms_kept": TERMS, "track_counts": track_counts,
                                   "sdk_version": sdk, "elapsed_seconds": elapsed})
        done += 1
        if done % REPORT_EVERY == 0:
            rate = done / max(time.time() - t_start, 1e-9)
            left = (len(regions) - skipped - done) / max(rate, 1e-9) / 3600
            la.log(f"{done} retrieved, {skipped} already archived, {i + 1}/{len(regions)} seen; "
                   f"{rate * 3600:.0f}/h, ~{left:.1f} h remaining")
            json.dump({"retrieved": done, "already_archived": skipped, "seen": i + 1, "total": len(regions),
                       "regions_per_hour": rate * 3600, "hours_remaining_estimate": left},
                      (TABLES / "saturation_retrieval_progress.json").open("w"), indent=1)
    json.dump({"retrieved": done, "already_archived": skipped, "total": len(regions), "complete": True},
              (TABLES / "saturation_retrieval_progress.json").open("w"), indent=1)
    la.log(f"saturation retrieval done: {done} retrieved, {skipped} already archived, {len(regions)} total")


if __name__ == "__main__":
    main()
