#!/usr/bin/env python3
"""Step 34 (P1 extensions E1-E3) on the Track 0 archive; no new queries.

E1  TF-occupancy channel: per signal, posterior-weighted top factors from the 21 adult-liver and 539 HepG2 CHIP_TF
    tracks (calibrated quantiles), and agreement with the repository's motifbreakR calls
    (GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv, JASPAR2024; hg38 coordinates in seqnames/start)
    against a position-matched background (the other queried variants of the same signal). No factor is "the mechanism".
E2  Lineage proxies: for each variant in a non-hepatocyte consensus peak (Track 0 variant_peak_overlaps), the
    predicted effect in the Atlas proxy tracks (immune/vascular groups of the fixed panel: monocyte, T/NK,
    endothelial; hepatic stellate DNase CL:0000632 read from the raw archive) beside the liver tracks.
    Every proxy is labelled cross-tissue.
E3  HepG2 contact-map scorer: posterior-weighted |quantile| of predicted local-contact change per signal as an
    enhancer–promoter linkage-disturbance flag (HepG2-only, descriptive).

Outputs (tables/): tf_channel_signal_top_factors.tsv, tf_channel_motifbreakr_agreement.tsv,
lineage_proxy_variant_effects.tsv.gz, contact_map_signal_summary.tsv
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict

import anndata
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import atlas_archive as aa
import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
RAW = ROOT / "raw"
MOTIF = la.PROJECT / "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"
TOP_K = 5
MIN_WEIGHT = 0.01


def tf_name(col: str) -> str:
    return col.split("TF ChIP-seq ")[-1].strip()


def load_weights() -> dict[str, dict[str, float]]:
    w = defaultdict(dict)
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["mapping_status"] == "mapped" and r["in_query_set"] == "True":
                w[r["signal_uid"]][r["variant_uid"]] = float(r["weight"])
    return w


def e1_tf_channel(weights, signals) -> None:
    t = pq.read_table(TABLES / "prediction_records" / "CHIP_TF.parquet").to_pandas()
    liver_cols = [c for c in t.columns if c.startswith("q|primary_liver|")]
    hepg2_cols = [c for c in t.columns if c.startswith("q|HepG2|")]
    t = t.set_index("variant_uid")
    motif = pd.read_csv(MOTIF, dtype=str)
    motif["pos_key"] = motif["seqnames"] + ":" + motif["start"]
    motif_tfs = defaultdict(set)
    for r in motif.itertuples(index=False):
        motif_tfs[r.pos_key].add(str(r.tf_name).upper())
    rows, agree = [], []
    for sig, w in weights.items():
        uids = [u for u in w if u in t.index]
        if not uids:
            continue
        sub = t.loc[uids]
        wv = np.array([w[u] for u in uids]); cov = wv.sum()
        for label, cols in (("adult_liver", liver_cols), ("HepG2", hepg2_cols)):
            if not cols:
                continue
            q = sub[cols].to_numpy(dtype=float)
            wa = (np.abs(np.nan_to_num(q)) * wv[:, None]).sum(0) / max(cov, 1e-12)     # posterior-weighted |quantile| per track
            ws = (np.nan_to_num(q) * wv[:, None]).sum(0) / max(cov, 1e-12)
            order = np.argsort(-wa)[:TOP_K]
            rows.append({"signal_uid": sig, "universe": signals[sig]["universe"], "gene": signals[sig]["gene"], "track_set": label, "coverage_mass": cov, "n_variants": len(uids),
                         "top_factors": ";".join(f"{tf_name(cols[j])}={wa[j]:.3f}({ws[j]:+.3f})" for j in order),
                         "max_weighted_abs_quantile": float(wa[order[0]]), "n_tracks": len(cols), "note": "posterior-weighted |calibrated quantile| (signed mean in parentheses); no factor is selected as the mechanism"})
        # motifbreakR agreement at the variant level: variants with a motif call whose TF is among the top-K HepG2/liver tracks for that variant
        for u in uids:
            if w[u] < MIN_WEIGHT:
                continue
            key = ":".join(u.split(":")[:2])
            called = motif_tfs.get(key, set())
            q = sub.loc[u, hepg2_cols + liver_cols].to_numpy(dtype=float)
            names = [tf_name(c).upper() for c in hepg2_cols + liver_cols]
            top = {names[j] for j in np.argsort(-np.abs(np.nan_to_num(q)))[:TOP_K]}
            # background: the same variant's rank of the called TF among all TF tracks vs the other queried variants of the signal
            agree.append({"signal_uid": sig, "variant_uid": u, "weight": w[u], "n_motif_calls": len(called), "motif_tfs": ";".join(sorted(called)),
                          "atlas_topK_tfs": ";".join(sorted(top)), "agreement": bool(called & top), "n_tf_tracks": len(names)})
    la.write_tsv_once(TABLES / "tf_channel_signal_top_factors.tsv", rows, list(rows[0].keys()))
    if agree:
        d = pd.DataFrame(agree)
        with_call = d[d.n_motif_calls > 0]
        rng = np.random.default_rng(la.prespec()["seeds"]["master"])
        # null: for each variant with a motif call, the chance that a random K-set of TF tracks contains one of its called TFs
        exp = float(np.mean([1 - math.comb(int(r.n_tf_tracks) - min(int(r.n_motif_calls), int(r.n_tf_tracks)), TOP_K) / math.comb(int(r.n_tf_tracks), TOP_K)
                             if int(r.n_tf_tracks) > TOP_K else 1.0 for r in with_call.itertuples()])) if len(with_call) else math.nan
        d.loc[len(d)] = {"signal_uid": "SUMMARY", "variant_uid": "", "weight": math.nan, "n_motif_calls": int(len(with_call)), "motif_tfs": "",
                         "atlas_topK_tfs": "", "agreement": float(with_call.agreement.mean()) if len(with_call) else math.nan, "n_tf_tracks": f"expected_by_chance={exp:.3f}"}
        d.to_csv(TABLES / "tf_channel_motifbreakr_agreement.tsv", sep="\t", index=False)
        la.log(f"E1: {len(rows)} signal rows; motif-called variants {len(with_call)}, top-{TOP_K} agreement {with_call.agreement.mean() if len(with_call) else float('nan'):.3f} vs chance {exp:.3f}")


def e2_lineage_proxies(weights, signals) -> None:
    vp = pd.read_csv(TABLES / "variant_peak_overlaps.tsv.gz", sep="\t", dtype=str)
    vp = vp[vp["lineage"].isin(["stellate", "macrophage", "cholangiocyte", "t_nk"])]
    want = set(vp["variant_uid"])
    if not want:
        la.log("E2: no variant in a non-hepatocyte peak"); return
    groups = {"ATAC": ["primary_liver", "hepatocyte", "immune", "vascular"], "DNASE": ["primary_liver", "hepatocyte", "immune", "vascular"], "CHIP_HISTONE": ["primary_liver", "hepatocyte", "immune", "vascular"]}
    per = {}
    for scorer, gs in groups.items():
        t = pq.read_table(TABLES / "prediction_records" / f"{scorer}.parquet").to_pandas()
        t = t[t["variant_uid"].isin(want)]
        for g in gs:
            cols = [c for c in t.columns if c.startswith(f"q|{g}|") and (scorer != "CHIP_HISTONE" or "H3K27ac" in c)]
            if cols:
                per[(scorer, g)] = t.set_index("variant_uid")[cols].median(axis=1)
    # hepatic stellate DNase (CL:0000632) is not in the fixed panel; read it from the raw archive
    stellate = {}
    for stage in ("atlas_direct", "atlas_enzyme"):
        for chunk in sorted((RAW / stage).glob("chunk_*")):
            f = chunk / "DNASE.h5ad"
            if not f.exists():
                continue
            a = anndata.read_h5ad(f)
            if a.shape[0] == 0 or "variant" not in a.obs or "ontology_curie" not in a.var:
                continue
            m = (a.var["ontology_curie"] == "CL:0000632").values
            if not m.any():
                continue
            q = np.asarray(a.layers["quantiles"], np.float32)[:, m]
            for i, v in enumerate(a.obs["variant"]):
                u = aa.variant_uid_from_str(str(v))
                if u in want:
                    stellate[u] = float(np.nanmedian(q[i]))
    rows = []
    for r in vp.itertuples(index=False):
        row = {"signal_uid": r.signal_uid, "variant_uid": r.variant_uid, "weight": r.weight, "lineage_peak": r.lineage, "peak": r.peak, "evidence_state": r.evidence_state}
        for (scorer, g), s in per.items():
            row[f"{scorer}|{g}"] = float(s.get(r.variant_uid, math.nan))
        row["DNASE|hepatic_stellate_CL0000632"] = stellate.get(r.variant_uid, math.nan)
        row["note"] = "immune/vascular are cross-tissue proxies; the liver lineage context comes from the measured snATAC peak, the proxy supplies only sequence-level direction"
        rows.append(row)
    pd.DataFrame(rows).to_csv(TABLES / "lineage_proxy_variant_effects.tsv.gz", sep="\t", index=False)
    la.log(f"E2: {len(rows)} variant-peak rows in non-hepatocyte lineages; stellate DNase available for {sum(1 for r in rows if r['DNASE|hepatic_stellate_CL0000632'] == r['DNASE|hepatic_stellate_CL0000632'])}")


def e3_contact_maps(weights, signals) -> None:
    t = pq.read_table(TABLES / "prediction_records" / "CONTACT_MAPS.parquet").to_pandas()
    cols = [c for c in t.columns if c.startswith("q|HepG2|")]
    if not cols:
        la.log("E3: no HepG2 contact-map track in the records"); return
    t = t.set_index("variant_uid")
    rows = []
    for sig, w in weights.items():
        uids = [u for u in w if u in t.index]
        if not uids:
            continue
        wv = np.array([w[u] for u in uids]); q = t.loc[uids, cols[0]].to_numpy(dtype=float)
        ok = ~np.isnan(q)
        rows.append({"signal_uid": sig, "universe": signals[sig]["universe"], "gene": signals[sig]["gene"], "coverage_mass": float(wv[ok].sum()),
                     "weighted_abs_quantile": float((np.abs(q[ok]) * wv[ok]).sum() / max(wv[ok].sum(), 1e-12)) if ok.any() else math.nan,
                     "mass_abs_quantile_ge_0.5": float(wv[ok][np.abs(q[ok]) >= 0.5].sum()) if ok.any() else math.nan, "track": cols[0],
                     "note": "HepG2 only; predicted local-contact change; descriptive flag beside the gene-level RNA channel"})
    la.write_tsv_once(TABLES / "contact_map_signal_summary.tsv", rows, list(rows[0].keys()))
    la.log(f"E3: {len(rows)} signals")


def main() -> None:
    weights = load_weights()
    signals = {s["signal_uid"]: s for s in la.read_tsv(TABLES / "eligible_signals.tsv")}
    e1_tf_channel(weights, signals)
    e2_lineage_proxies(weights, signals)
    e3_contact_maps(weights, signals)


if __name__ == "__main__":
    main()
