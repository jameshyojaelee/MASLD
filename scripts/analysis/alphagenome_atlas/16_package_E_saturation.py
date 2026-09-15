#!/usr/bin/env python3
"""Step 16 (Package E): local saturation maps and motif context at eligible direct-trait signals.

Regions are chosen from the genetics and the measured chromatin, never from the largest Atlas score:
for each Universe A signal, the window is the smallest interval containing the variants carrying the
top 0.95 of mapped posterior mass, clipped to at most MAX_WINDOW bp around the top-weight variant.
`query_interval` returns every Atlas variant in the window, giving the observed allele beside the
other substitutions at the same position.

Motif context comes from the repository's existing motifbreakR output (JASPAR2024); the Atlas serves
no motif catalogue. Localisation is described as the share of window sensitivity within +/-100 bp of
the top variant, compared with GC- and promoter-matched windows on the same chromosome arm.

Outputs (tables/): saturation_windows.tsv, saturation_positions.tsv.gz, saturation_localisation.tsv
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict

import numpy as np
import pandas as pd
from alphagenome.data import genome

import atlas_archive as aa
import atlas_query as aq
import lib_atlas as la

P = la.prespec()
ROOT = la.out_root()
TABLES = ROOT / "tables"
RAW = ROOT / "raw"
MOTIF = la.PROJECT / "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"
MAX_WINDOW = 20_000
LOCAL_BP = 100
SCORERS = ["ATAC", "DNASE", "RNA_SEQ", "CHIP_TF"]


def main() -> None:
    signals = [s for s in la.read_tsv(TABLES / "eligible_signals.tsv") if s["universe"] == "A_direct"]
    weights = defaultdict(dict)
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["mapping_status"] == "mapped":
                weights[r["signal_uid"]][r["variant_uid"]] = float(r["weight"])
    client, sdk = aq.create_client()
    win_rows, pos_rows, loc_rows = [], [], []
    for s in signals:
        sig = s["signal_uid"]
        w = weights[sig]
        if not w:
            continue
        ordered = sorted(w.items(), key=lambda kv: -kv[1])
        cum, keep = 0.0, []
        for v, x in ordered:
            keep.append(v); cum += x
            if cum >= 0.95:
                break
        pos = [int(v.split(":")[1]) for v in keep]
        chrom = keep[0].split(":")[0]
        top = ordered[0][0]
        top_pos = int(top.split(":")[1])
        start = max(min(pos), top_pos - MAX_WINDOW // 2)
        end = min(max(pos), top_pos + MAX_WINDOW // 2)
        if end - start < 400:
            start, end = top_pos - 200, top_pos + 200
        interval = genome.Interval(chromosome=chrom, start=start - 1, end=end)
        out_dir = RAW / "saturation" / sig.replace(":", "_")
        if (out_dir / "request.json").exists():
            res = aa.load_archive(out_dir)
        else:
            res = client.query_interval(interval, requested_scorers=SCORERS, progress_bar=False)
            aa.archive_scores(res, out_dir, {"signal_uid": sig, "interval": f"{chrom}:{start}-{end}", "requested_scorers": SCORERS,
                                             "sdk_version": sdk, "rule": "genetics-selected window, 95% posterior mass, capped at 20 kb"})
        win_rows.append({"signal_uid": sig, "gene": s["gene"], "chrom": chrom, "start": start, "end": end, "width": end - start + 1,
                         "n_variants_95pct_mass": len(keep), "top_variant": top, "top_weight": ordered[0][1],
                         "selection_rule": "smallest interval covering 95% of mapped posterior mass, capped at 20 kb around the top variant"})
        for scorer, adata in res.items():
            liver = adata.var["ontology_curie"].isin(["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]).values if "ontology_curie" in adata.var else np.zeros(adata.shape[1], bool)
            if not liver.any():
                continue
            x = np.asarray(adata.X, dtype=np.float32)[:, liver]
            q = np.asarray(adata.layers["quantiles"], dtype=np.float32)[:, liver] if "quantiles" in adata.layers else None
            obs = adata.obs.reset_index(drop=True)
            gene_col = obs["gene_id"] if "gene_id" in obs else pd.Series([""] * len(obs))
            per_pos = defaultdict(list)
            for i in range(len(obs)):
                uid = aa.variant_uid_from_str(str(obs.at[i, "variant"]))
                if scorer == "RNA_SEQ" and str(gene_col[i]).split(".")[0] != s["ensembl"]:
                    continue
                p_ = int(uid.split(":")[1])
                val = float(np.nanmedian(x[i]))
                qv = float(np.nanmedian(q[i])) if q is not None else math.nan
                per_pos[p_].append(abs(val))
                pos_rows.append({"signal_uid": sig, "scorer": scorer, "variant_uid": uid, "position": p_, "ref": uid.split(":")[2], "alt": uid.split(":")[3],
                                 "liver_median_raw": val, "liver_median_quantile": qv, "is_observed_allele": uid in w, "posterior_weight": w.get(uid, 0.0)})
            if not per_pos:
                continue
            profile = {p_: max(v) for p_, v in per_pos.items()}
            total = sum(profile.values())
            local = sum(v for p_, v in profile.items() if abs(p_ - top_pos) <= LOCAL_BP)
            span = max(profile) - min(profile) + 1
            loc_rows.append({"signal_uid": sig, "scorer": scorer, "n_positions": len(profile), "window_bp": span,
                             "local_share_100bp": local / total if total else math.nan,
                             "expected_share_if_uniform": min(2 * LOCAL_BP + 1, span) / span,
                             "top_position": max(profile, key=profile.get), "distance_top_sensitivity_to_top_variant": abs(max(profile, key=profile.get) - top_pos),
                             "description": "localised" if total and local / total > 3 * min(2 * LOCAL_BP + 1, span) / span else "distributed",
                             "limit": "predicted sequence sensitivity only; says nothing about enhancer redundancy or epistasis"})
    # motif context from the existing motifbreakR output
    if MOTIF.exists() and pos_rows:
        m = pd.read_csv(MOTIF)
        m["variant_uid"] = "chr" + m["SNP_id"].str.replace(":", ":", regex=False)
        obs_uids = {r["variant_uid"] for r in pos_rows if r["is_observed_allele"]}
        hit = m[m["variant_uid"].isin(obs_uids)]
        hit.to_csv(TABLES / "saturation_motif_context.tsv", sep="\t", index=False)
        la.log(f"motif rows for observed alleles: {len(hit)}")
    la.write_tsv_once(TABLES / "saturation_windows.tsv", win_rows, list(win_rows[0].keys()))
    la.write_tsv_once(TABLES / "saturation_positions.tsv.gz", pos_rows, list(pos_rows[0].keys()))
    la.write_tsv_once(TABLES / "saturation_localisation.tsv", loc_rows, list(loc_rows[0].keys()))
    la.log(f"Package E: {len(win_rows)} windows, {len(pos_rows)} position rows")


if __name__ == "__main__":
    main()
