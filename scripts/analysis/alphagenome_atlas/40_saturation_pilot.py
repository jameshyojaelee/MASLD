#!/usr/bin/env python3
"""Step 40 (P4 pilot): time, payload and archive size of `query_interval` saturation on measured peaks.

Retrieval is NOT ontology-filtered. Passing the 50-curie liver+proxy list to query_interval makes the call
fail with UNAVAILABLE after 60-208 s, while the identical request without it returns in 2.2 s; the ontology
filter is therefore applied client-side at analysis time, as every other step already does.

Fixed before any result: the region universes (H3K27ac 96,460 regions of the released predictor; snATAC consensus
peaks restricted to DA / genetics / program-linked peaks in step 41), the ontology filter (liver + HepG2 + the four
cross-tissue proxies), and the scorers. The pilot takes the first N_PILOT regions of each universe in genomic order,
records seconds and bytes per region, whether interval rows include indels, and the number of variants per region.

Outputs: raw/saturation_pilot/<region>/, tables/saturation_pilot_report.json
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np
from alphagenome.data import genome

import atlas_archive as aa
import atlas_query as aq
import lib_atlas as la

ROOT = la.out_root()
RAW = ROOT / "raw"
TABLES = ROOT / "tables"
BENCH = la.PROJECT / "Analysis/MASLD_Model_Benchmark"
RELEASE = BENCH / "release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz"
PEAKS = la.PROJECT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/consensus_peak_manifest.tsv"
N_PILOT = int(sys.argv[1]) if len(sys.argv) > 1 else 20
SCORERS = ["ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE", "AVI_SCORE"]
ONTOLOGY = {
    "liver": ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182", "CL:0000632", "CL:0002571"],
    "HepG2": ["EFO:0001187"],
    "fibroblast_proxy": ["CL:0002547", "CL:0002548", "CL:0002549", "CL:0002550", "CL:0002551", "CL:0002552", "CL:0002553", "CL:0002554", "CL:0002555", "CL:0002556", "CL:0002557", "CL:0002558"],
    "monocyte_macrophage_proxy": ["CL:0000235", "CL:0000576", "CL:0000862", "CL:0000863", "CL:0001054", "CL:0002057", "CL:0002397"],
    "endothelial_proxy": ["CL:0000115", "CL:0002138", "CL:0002144", "CL:0002188", "CL:0002543", "CL:0002544", "CL:0002545", "CL:0002618", "CL:1000413", "CL:1000892", "CL:1001568", "CL:2000010"],
    "t_nk_proxy": ["CL:0000084", "CL:0000097", "CL:0000136", "CL:0000351", "CL:0000623", "CL:0000624", "CL:0000625", "CL:0000792", "CL:0000798", "CL:0000823", "CL:0000895", "CL:0000897"],
}
TERMS = sorted({c for v in ONTOLOGY.values() for c in v})


def universes() -> dict[str, list[tuple[str, int, int]]]:
    rel = np.load(RELEASE, allow_pickle=True)
    keys = [str(k) for k in rel["prof_region_key"]]
    h3 = []
    for k in keys:
        chrom, rng_ = k.split(":"); a, b = rng_.split("-")
        h3.append((chrom, int(a), int(b)))
    peaks = []
    for r in la.read_tsv(PEAKS):
        if r["lineage"] == "hepatocyte" and r["blacklist_overlap"] != "TRUE" and r["supported_both"] == "TRUE":
            peaks.append((r["chrom"], int(r["start0"]), int(r["end"])))
    return {"h3k27ac_regions": h3[:N_PILOT], "snatac_hepatocyte_peaks_supported_both": peaks[:N_PILOT]}


def main() -> None:
    client, sdk = aq.create_client()
    report = {"ontology_terms_applied": "client-side after retrieval, not server-side (see module docstring)",
              "ontology_terms": TERMS, "scorers": SCORERS, "sdk_version": sdk, "regions": []}
    for uni, regs in universes().items():
        for chrom, start0, end in regs:
            key = f"{chrom}:{start0}-{end}"
            d = RAW / "saturation_pilot" / uni / key.replace(":", "_")
            interval = genome.Interval(chromosome=chrom, start=start0, end=end)
            t0 = time.time()
            if (d / "request.json").exists():
                res = aa.load_archive(d); elapsed = json.load((d / "request.json").open()).get("elapsed_seconds")
            else:
                # ontology_terms is NOT passed: a 50-curie filter makes query_interval fail with UNAVAILABLE
                # after 60-208 s, while the identical request without it returns in 2.2 s. The ontology
                # filter is applied client-side after retrieval instead, which is where every other step
                # already applies it (liver_mask / liver_track_mask). Retrieval keeps all tracks.
                res = aq.call_with_quota_retry(lambda: client.query_interval(interval, requested_scorers=SCORERS, progress_bar=False, max_workers=4), label=key)
                elapsed = time.time() - t0
                res, track_counts = aq.filter_tracks_to_ontology(res, TERMS)
                aa.archive_scores(res, d, {"universe": uni, "interval": key, "requested_scorers": SCORERS,
                                           "ontology_terms_requested": None, "ontology_filter": "archive-time, client-side",
                                           "ontology_terms_kept": TERMS, "track_counts": track_counts,
                                           "sdk_version": sdk, "elapsed_seconds": elapsed})
            size = sum(f.stat().st_size for f in d.glob("*"))
            entry = {"universe": uni, "region": key, "width": end - start0, "seconds": elapsed, "archive_bytes": size}
            for scorer, a in res.items():
                n_var = int(a.shape[0]); n_tr = int(a.shape[1])
                indel = 0
                if n_var and "variant" in a.obs:
                    vs = [str(v) for v in a.obs["variant"]]
                    indel = sum(1 for v in vs if len(v.split(":")[-1].split(">")[0]) != 1 or len(v.split(">")[-1]) != 1)
                entry[scorer] = {"n_rows": n_var, "n_tracks": n_tr, "n_indel_rows": indel}
            report["regions"].append(entry)
            la.log(f"{uni} {key}: {elapsed:.1f}s, {size/1e6:.2f} MB, ATAC rows {entry.get('ATAC', {}).get('n_rows')}")
    secs = [r["seconds"] for r in report["regions"] if r["seconds"]]
    mbs = [r["archive_bytes"] / 1e6 for r in report["regions"]]
    report["summary"] = {"n_regions": len(report["regions"]), "median_seconds": float(np.median(secs)) if secs else None, "median_MB": float(np.median(mbs)) if mbs else None,
                         "projected_hours_96460": (np.median(secs) * 96460 / 3600) if secs else None, "projected_GB_96460": (np.median(mbs) * 96460 / 1000) if mbs else None,
                         "any_indel_rows": any(v.get("n_indel_rows", 0) for r in report["regions"] for v in r.values() if isinstance(v, dict))}
    json.dump(report, (TABLES / "saturation_pilot_report.json").open("w"), indent=1)
    la.log(f"pilot summary: {report['summary']}")


if __name__ == "__main__":
    main()
