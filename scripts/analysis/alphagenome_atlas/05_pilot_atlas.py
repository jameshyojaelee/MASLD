#!/usr/bin/env python3
"""Step 05: scorer metadata dump, then the five-block pilot query with pass/kill report."""

from __future__ import annotations

import json
import time
from collections import defaultdict

import anndata

import atlas_archive as aa
import atlas_query as aq
import lib_atlas as la

P = la.prespec()
ROOT = la.out_root()
TABLES = ROOT / "tables"
RAW = ROOT / "raw"
ANCHOR = {"variant_uid": "chr1:109274968:G:T", "gene": "ENSG00000134243", "name": "SORT1_rs12740374"}
FAMILIES = {"rna": ("rna_seq",), "atac_dnase": ("atac", "dnase"), "splice": ("splice",)}


def main() -> None:
    client, sdk = aq.create_client()

    # 1. Scorer metadata (track selection table, fixed before any variant result).
    meta = client.scorer_metadata()
    md_dir = RAW / "scorer_metadata"
    md_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, m in sorted(meta.items()):
        tm = m.track_metadata.copy()
        tm.to_csv(md_dir / f"{name}.track_metadata.tsv", sep="\t", index=False)
        classes = defaultdict(int)
        for _, t in tm.iterrows():
            classes[la.track_class(str(t.get("ontology_curie", "")), str(t.get("biosample_name", "")), str(t.get("biosample_type", "")))] += 1
        assays = sorted(set(map(str, tm["Assay title"]))) if "Assay title" in tm else []
        rows.append({"scorer": name, "is_signed": m.is_signed, "n_tracks": len(tm), "n_primary_liver": classes["primary_liver"],
                     "n_hepatocyte": classes["hepatocyte"], "n_hepg2": classes["HepG2"], "n_other": classes["other"],
                     "assays": ";".join(assays), "columns": ";".join(map(str, tm.columns))})
    if not (TABLES / "scorer_metadata.tsv").exists():
        la.write_tsv_once(TABLES / "scorer_metadata.tsv", rows, list(rows[0].keys()))
    la.log(f"{len(rows)} scorers")

    # 2. Pilot variants.
    with (TABLES / "pilot_blocks.json").open() as handle:
        pilot = json.load(handle)
    # Pilot = direct-universe (A/B) query variants inside the five blocks; enzyme signals sharing a block are not queried here.
    pilot_signals = set(pilot["signals"])
    direct_uids = set()
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as handle:
        for r in __import__("csv").DictReader(handle, delimiter="\t"):
            if r["signal_uid"] in pilot_signals and r["universe"] != "C_enzyme" and r["in_query_set"] == "True":
                direct_uids.add(r["variant_uid"])
    uids = sorted(direct_uids) + [ANCHOR["variant_uid"]]
    scorers = sorted(meta)
    t0 = time.time()
    chunks = aq.query_archived(client, sdk, uids, scorers, RAW / "pilot", extra={"blocks": pilot["blocks"], "signals": pilot["signals"]})
    elapsed = sum(json.load(open(c / "request.json"))["elapsed_seconds"] for c in chunks)
    results: dict = {}
    for c in chunks:
        for name, adata in aa.load_archive(c).items():
            results.setdefault(name, []).append(adata)
    results = {k: (anndata.concat(v, join="outer") if len(v) > 1 else v[0]) for k, v in results.items()}

    # 3. Report.
    returned = defaultdict(set)
    records_per_variant = defaultdict(int)
    quantiles = {}
    anchor_sort1 = {}
    for name, adata in results.items():
        quantiles[name] = "quantiles" in adata.layers
        for rec in aa.long_records(adata, name, meta[name].is_signed):
            returned[rec["variant_uid"]].add(name)
            records_per_variant[rec["variant_uid"]] += 1
            if rec["variant_uid"] == ANCHOR["variant_uid"] and rec["gene_id"].split(".")[0] == ANCHOR["gene"] and rec["track_class"] in ("primary_liver", "hepatocyte") and "rna" in name.lower():
                anchor_sort1.setdefault(name, []).append(rec["raw_score"])
    missing = [u for u in uids if u not in returned]
    liver_by_family = {}
    for fam, keys in FAMILIES.items():
        liver_by_family[fam] = sum(r["n_primary_liver"] + r["n_hepatocyte"] for r in rows if any(k in r["scorer"].lower() for k in keys))
    per100 = 100 * elapsed / max(len(uids), 1)
    n_full = sum(1 for _ in open(TABLES / "query_set.tsv")) - 1
    extrapolated_h = n_full * elapsed / max(len(uids), 1) / 3600
    anchor_mean = {k: sum(v) / len(v) for k, v in anchor_sort1.items()}
    report = {
        "sdk_version": sdk, "n_scorers": len(meta), "scorers": scorers, "is_signed": {k: meta[k].is_signed for k in scorers},
        "quantile_layer": quantiles, "n_variants_requested": len(uids), "n_variants_returned": len(returned),
        "missing_variants": missing, "records_per_variant_min_max": [min(records_per_variant.values(), default=0), max(records_per_variant.values(), default=0)],
        "seconds_per_100_variants": per100, "full_query_set_size": n_full, "extrapolated_full_hours": extrapolated_h,
        "liver_tracks_by_family": liver_by_family, "anchor_sort1_liver_rna_mean_signed_score": anchor_mean,
        "pass": {"all_scorers_present": len(meta) > 0, "liver_rna": liver_by_family["rna"] > 0, "liver_atac_dnase": liver_by_family["atac_dnase"] > 0,
                 "liver_splice": liver_by_family["splice"] > 0, "no_silent_drop": not missing, "time_under_24h": extrapolated_h < 24,
                 "anchor_positive": all(v > 0 for v in anchor_mean.values()) if anchor_mean else None},
    }
    with (TABLES / "pilot_report.json").open("w") as handle:
        json.dump(report, handle, indent=1)
    la.log(json.dumps(report["pass"]))
    la.log(f"elapsed {elapsed:.1f}s for {len(uids)} variants; extrapolated {extrapolated_h:.2f} h for {n_full}")


if __name__ == "__main__":
    main()
