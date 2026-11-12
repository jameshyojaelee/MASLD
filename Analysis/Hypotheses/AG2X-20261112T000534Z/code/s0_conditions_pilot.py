#!/usr/bin/env python3
"""Stage 0, steps 2-3: dump the human condition catalog, flag liver-relevant conditions, and run the
200-call pilot on A1 caQTL leads in the bulk-liver condition (20 rescored to test determinism).

Bulk liver for the pilot is the condition whose ontology is exactly UBERON:0002107 (liver). If it
is absent the pilot is skipped and the catalog alone is written. Pilot uses seed 20261111 to draw
200 of the 3,845 tier-A4 leads.
"""
from __future__ import annotations

import json
import re
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

OUT = L.PKG / "stage0"
A1 = L.REPO / "GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z/tables/endpoint1_matched_tierA4.tsv.gz"
PILOT_SCORERS = ["ATAC", "DNASE", "SC_ATAC", "RNA_SEQ", "SC_GENE_EXPRESSION"]

TERMS = {
    "liver_bulk": r"UBERON:0002107|UBERON:0001114|UBERON:0001115|\bliver\b",
    "hepatocyte": r"CL:0000182|hepatocyte",
    "stellate": r"CL:0000632|stellate",
    "kupffer_macrophage": r"CL:0000091|kupffer|macrophage",
    "sinusoidal_endothelial": r"CL:1000398|CL:0002262|sinusoid",
    "cholangiocyte": r"CL:1000488|cholangiocyte|bile duct|biliary",
    "hepg2": r"EFO:0001187|hepg2",
    "cd14_monocyte": r"CL:0001054|CL:0002057|CD14",
    "neuron": r"CL:0000540|neuron",
    "skeletal_muscle": r"UBERON:0001134|UBERON:0011907|CL:0000187|CL:0000188|skeletal muscle|myocyte|myotube",
}


def main():
    t_all = time.time()
    counter = {}
    client = L.create_v2()
    from alphagenome.models.v2 import dna_client
    cond = client.supported_conditions(dna_client.Organism.HOMO_SAPIENS)
    cond.to_csv(OUT / "supported_conditions_human.tsv", sep="\t", index=False)
    L.log(f"supported_conditions: {len(cond)} rows, columns {list(cond.columns)}")
    try:
        meta = client.observed_metadata(dna_client.Organism.HOMO_SAPIENS)
        rows = []
        for ot in meta.__dict__ if hasattr(meta, "__dict__") else []:
            pass
        # OutputMetadata exposes one DataFrame per output type
        for name in dir(meta):
            if name.startswith("_"):
                continue
            v = getattr(meta, name)
            if isinstance(v, pd.DataFrame):
                d = v.copy()
                d.insert(0, "output_type", name)
                rows.append(d)
        if rows:
            om = pd.concat(rows, ignore_index=True)
            om.to_csv(OUT / "observed_metadata_human.tsv.gz", sep="\t", index=False)
            L.log(f"observed_metadata: {len(om)} track rows")
    except Exception as e:  # noqa: BLE001
        L.log(f"observed_metadata failed: {type(e).__name__}: {e}")

    txt = (cond["condition"].astype(str) + " | " + cond["biosample_name"].astype(str))
    flagged = []
    for cat, rx in TERMS.items():
        m = txt.str.contains(rx, flags=re.IGNORECASE, regex=True)
        d = cond[m].copy()
        d.insert(0, "category", cat)
        flagged.append(d)
    flagged = pd.concat(flagged, ignore_index=True)
    flagged.to_csv(OUT / "conditions_flagged.tsv", sep="\t", index=False)
    L.log(f"flagged {len(flagged)} condition rows across {flagged.category.nunique()} categories")

    bulk = cond[cond["condition"].astype(str) == "UBERON:0002107"]
    summary = {"n_conditions": int(len(cond)), "package_version": L.package_version(), "utc": L.utc(),
               "bulk_liver_present": bool(len(bulk))}
    if not len(bulk):
        L.log("no exact UBERON:0002107 condition; pilot skipped")
        json.dump(summary, open(OUT / "stage0_summary.json", "w"), indent=1)
        return

    a1 = pd.read_csv(A1, sep="\t", usecols=["key", "chr_x", "pos_hg38", "ref", "alt"])
    a1["uid"] = [f"{c}:{p}:{r}:{a}" for c, p, r, a in zip(a1.chr_x, a1.pos_hg38, a1.ref, a1.alt)]
    rng = np.random.default_rng(L.SEED)
    pick = sorted(rng.choice(len(a1), size=200, replace=False).tolist())
    uids = a1.uid.iloc[pick].tolist()
    bad = L.assert_ref(uids)
    summary["pilot_ref_mismatches"] = len(bad)
    if bad:
        raise SystemExit(f"REF mismatch against FASTA: {bad[:5]}")
    pd.DataFrame({"uid": uids}).to_csv(OUT / "pilot_variants.tsv", sep="\t", index=False)

    # first call alone: if a scorer is rejected for this condition, find which one
    try:
        L.score_chunked(client, "v2", uids[:1], "UBERON:0002107", PILOT_SCORERS, OUT / "pilot_first", chunk=1,
                        max_workers=1, counter=counter)
        scorers = PILOT_SCORERS
    except Exception as e:  # noqa: BLE001
        L.log(f"first pilot call failed with all scorers: {type(e).__name__}: {getattr(e, 'details', lambda: e)()}")
        scorers = []
        for s in PILOT_SCORERS:
            try:
                L.score_chunked(client, "v2", uids[:1], "UBERON:0002107", [s], OUT / f"pilot_probe_{s}", chunk=1,
                                max_workers=1, counter=counter)
                scorers.append(s)
            except Exception as e2:  # noqa: BLE001
                L.log(f"scorer {s} rejected: {type(e2).__name__}: {getattr(e2, 'details', lambda: e2)()}")
    summary["pilot_scorers_accepted"] = scorers

    t0 = time.time()
    L.score_chunked(client, "v2", uids[1:], "UBERON:0002107", scorers, OUT / "pilot_main", chunk=50,
                    max_workers=4, counter=counter)
    el = time.time() - t0
    summary["pilot_main_calls"] = len(uids) - 1
    summary["pilot_main_seconds"] = el
    summary["pilot_calls_per_second_workers4"] = (len(uids) - 1) / el

    rep = uids[1:21]
    L.score_chunked(client, "v2", rep, "UBERON:0002107", scorers, OUT / "pilot_repeat", chunk=20,
                    max_workers=4, counter=counter)
    a = pd.read_parquet(OUT / "pilot_main" / "chunk_00000.parquet")
    b = pd.read_parquet(OUT / "pilot_repeat" / "chunk_00000.parquet")
    keys = ["variant_id", "variant_scorer", "track_name", "track_strand", "gene_id", "junction_Start", "junction_End"]
    keys = [k for k in keys if k in a.columns and k in b.columns]
    a = a[a.variant_id.isin(b.variant_id.unique())]
    m = a.merge(b, on=keys, suffixes=("_a", "_b"))
    d = (m.raw_score_a - m.raw_score_b).abs()
    summary["determinism_rows_compared"] = int(len(m))
    summary["determinism_rows_first"] = int(len(a))
    summary["determinism_rows_repeat"] = int(len(b))
    summary["determinism_max_abs_diff"] = float(d.max()) if len(d) else None
    summary["determinism_median_abs_score"] = float(m.raw_score_a.abs().median()) if len(m) else None
    summary["determinism_sign_agreement"] = float((np.sign(m.raw_score_a) == np.sign(m.raw_score_b)).mean()) if len(m) else None
    sizes = [json.loads(x) for x in open(OUT / "pilot_main" / "chunk_log.jsonl")]
    summary["pilot_rows_per_variant"] = sum(s["rows"] for s in sizes) / sum(s["n_variants"] for s in sizes)
    summary["pilot_parquet_bytes_per_variant"] = sum(s["bytes_parquet"] for s in sizes) / sum(s["n_variants"] for s in sizes)
    main0 = pd.read_parquet(OUT / "pilot_main" / "chunk_00000.parquet")
    summary["pilot_rows_by_scorer_and_mode"] = (main0.groupby(["output_type", "imputation_mode"]).size()
                                                .rename("rows").reset_index().to_dict(orient="records"))
    summary["hosted_calls"] = counter
    summary["wall_seconds"] = time.time() - t_all
    json.dump(summary, open(OUT / "stage0_summary.json", "w"), indent=1, default=str)
    L.log(json.dumps(summary, default=str)[:3000])


if __name__ == "__main__":
    main()
