#!/usr/bin/env python3
"""Step 21 (P2): known-mechanism anchor panel — what the Atlas can and cannot see in liver.

Reads 21_anchor_panel_prespec.json (fixed before any score), resolves hg38 alleles (rsID -> hg19 from the SuSiE
aggregate -> repository liftOver rule -> FASTA-validated reference), point-queries every served anchor with all
scorers, retrieves the 2-kb window with query_interval as the position-matched background, and scores each anchor
against its prespecified expected channel.

Inputs : tables/anchor_hg38.tsv (from 21a_anchor_liftover.R), the prespec, the FASTA
Outputs: raw/atlas_anchor_points/, raw/atlas_anchor_windows/<rsid>/, tables/anchor_panel.tsv, tables/anchor_panel_summary.json
"""

from __future__ import annotations

import json
import math
import time
from collections import defaultdict

import anndata
import numpy as np
import pandas as pd
import pysam
from alphagenome.data import genome

import atlas_archive as aa
import atlas_query as aq
import lib_atlas as la

ROOT = la.out_root()
RAW = ROOT / "raw"
TABLES = ROOT / "tables"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
PRESPEC = json.load((la.SCRIPT_DIR / "21_anchor_panel_prespec.json").open())
LIVER = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
SCORERS = ["RNA_SEQ", "SPLICE_SITE_USAGE", "ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "AVI_SCORE", "CAGE"]
CHANNELS = {"rna_gene": "RNA_SEQ", "splice_site_usage": "SPLICE_SITE_USAGE", "atac": "ATAC", "dnase": "DNASE", "h3k27ac": "CHIP_HISTONE"}


# ---------------------------------------------------------------- helpers (tested)
def resolve_alleles(fasta_base: str, alleles: str) -> tuple[str, str] | None:
    """ref = the FASTA base; alt = the other reported allele. Returns None when the FASTA base is not among the alleles."""
    parts = [a.upper() for a in alleles.split("/")]
    base = fasta_base.upper()
    if base not in parts:
        return None
    others = [a for a in parts if a != base]
    if len(others) != 1:
        return None
    return base, others[0]


def rank_percentile(value: float, background: list[float]) -> float:
    """Share of background |values| that are <= |value| (1.0 = the anchor is the largest in its window)."""
    bg = [abs(b) for b in background if b == b]
    if not bg or value != value:
        return math.nan
    return sum(1 for b in bg if b <= abs(value)) / len(bg)


def score_expectation(expected_channel: str, requested: bool, returned: bool,
                      channel_values: dict, channel_ranks: dict) -> str:
    """Verdict for one anchor against its prespecified channel.

    `requested` is whether a queryable hg38 uid was built (an indel or an unresolved coordinate is not).
    `returned` is whether the Atlas actually gave back a row for that uid. The two are different: the GNMT
    haplotype pair resolved cleanly and was queried, and the Atlas returned nothing for either variant.
    Calling that a `miss` would report a negative about a variant that was never measured, so it is
    `uncovered`.
    """
    if not requested:
        return "not_served"
    if expected_channel == "none":
        return "not_applicable"
    if expected_channel in ("none_stated",):
        return "not_prespecified"
    if not returned or not any(v == v for v in channel_values.values()):
        return "uncovered"
    candidates = list(CHANNELS) if expected_channel == "any_regulatory" else [expected_channel]
    if expected_channel == "any_regulatory":
        candidates = [c for c in candidates if c != "rna_gene"] + ["rna_gene"]
    for c in candidates:
        v, r = channel_values.get(c, math.nan), channel_ranks.get(c, math.nan)
        if v == v and r == r and abs(v) >= 0.5 and r >= 0.9:
            return "hit"
    return "miss"


def liver_mask(var: pd.DataFrame, histone: str | None = None) -> np.ndarray:
    m = var["ontology_curie"].isin(LIVER).values if "ontology_curie" in var else np.zeros(len(var), bool)
    if histone is not None and "name" in var:
        m &= var["name"].astype(str).str.contains(histone).values
    return m


def channel_summary(res: dict, uid: str, ensembl: str) -> tuple[dict, dict]:
    """Per channel: liver median quantile at this uid (gene-scoped for RNA/splice); plus AVI and top TF tracks."""
    vals, extra = {}, {}
    for ch, scorer in CHANNELS.items():
        a = res.get(scorer)
        if a is None or a.shape[0] == 0 or "variant" not in a.obs:
            vals[ch] = math.nan
            continue
        obs = a.obs.reset_index(drop=True)
        rows = [i for i in range(len(obs)) if aa.variant_uid_from_str(str(obs.at[i, "variant"])) == uid
                and (scorer not in ("RNA_SEQ", "SPLICE_SITE_USAGE") or str(obs.at[i, "gene_id"]).split(".")[0] == ensembl)]
        if not rows:
            vals[ch] = math.nan
            continue
        mask = liver_mask(a.var, "H3K27ac" if ch == "h3k27ac" else None)
        if not mask.any():
            vals[ch] = math.nan
            continue
        q = np.asarray(a.layers["quantiles"], np.float32)[rows][:, mask]
        vals[ch] = float(np.nanmax(np.abs(q)) * np.sign(q.flat[int(np.nanargmax(np.abs(q)))])) if ch == "splice_site_usage" else float(np.nanmedian(q))
    a = res.get("AVI_SCORE")
    if a is not None and a.shape[0] and "variant" in a.obs:
        obs = a.obs.reset_index(drop=True)
        rows = [i for i in range(len(obs)) if aa.variant_uid_from_str(str(obs.at[i, "variant"])) == uid]
        if rows:
            extra["avi_raw"] = float(np.asarray(a.X)[rows[0], 0]); extra["avi_quantile"] = float(np.asarray(a.layers["quantiles"])[rows[0], 0]) if "quantiles" in a.layers else math.nan
    a = res.get("CHIP_TF")
    if a is not None and a.shape[0] and "variant" in a.obs:
        obs = a.obs.reset_index(drop=True)
        rows = [i for i in range(len(obs)) if aa.variant_uid_from_str(str(obs.at[i, "variant"])) == uid]
        if rows:
            var = a.var.reset_index(drop=True)
            keep = var["biosample_name"].astype(str).str.contains("liver|HepG2", case=False, regex=True).values
            q = np.asarray(a.layers["quantiles"], np.float32)[rows[0]]
            order = [j for j in np.argsort(-np.abs(np.where(keep, q, 0)))[:5] if keep[j]]
            extra["tf_top"] = ";".join(f"{var.at[j, 'name']}={q[j]:+.2f}" for j in order)
    return vals, extra


def main() -> None:
    fasta = pysam.FastaFile(FASTA)
    hg38 = {r["rsid"]: r for r in la.read_tsv(TABLES / "anchor_hg38.tsv")}
    client, sdk = aq.create_client()
    anchors = []
    for a in PRESPEC["anchors"]:
        row = {**a, "served": False, "variant_uid": "", "hg38_chrom": "", "hg38_pos": "", "ref": "", "alt": ""}
        h = hg38.get(a["rsid"])
        if h is None or int(h["n_map"]) != 1:
            row["allele_state"] = "no_unique_hg38_mapping"
            anchors.append(row); continue
        chrom, pos = h["hg38_chrom"], int(h["hg38_pos"])
        alleles = h["alleles"]
        if any(len(x) != 1 for x in alleles.split("/")):
            row.update({"hg38_chrom": chrom, "hg38_pos": pos, "allele_state": "indel_not_served_by_point_endpoint"})
            anchors.append(row); continue
        base = fasta.fetch(chrom, pos - 1, pos)
        ra = resolve_alleles(base, alleles)
        if ra is None:
            row.update({"hg38_chrom": chrom, "hg38_pos": pos, "allele_state": f"fasta_base_{base}_not_in_{alleles}"})
            anchors.append(row); continue
        ref, alt = ra
        row.update({"hg38_chrom": chrom, "hg38_pos": pos, "ref": ref, "alt": alt, "variant_uid": la.variant_uid(chrom, pos, ref, alt),
                    "allele_state": "palindromic_ref_from_fasta" if h.get("palindromic") == "TRUE" else "resolved", "served": True})
        anchors.append(row)
    uids = sorted({r["variant_uid"] for r in anchors if r["served"]})
    aq.query_archived(client, sdk, uids, SCORERS, RAW / "atlas_anchor_points", extra={"stage": "anchor_points"})
    point = {}
    for c in sorted((RAW / "atlas_anchor_points").glob("chunk_*")):
        point[c] = aa.load_archive(c)
    win = PRESPEC["window_bp"]
    out_rows = []
    for r in anchors:
        vals, ranks, extra = {c: math.nan for c in CHANNELS}, {c: math.nan for c in CHANNELS}, {}
        if r["hg38_chrom"]:
            wdir = RAW / "atlas_anchor_windows" / r["rsid"]
            interval = genome.Interval(chromosome=r["hg38_chrom"], start=int(r["hg38_pos"]) - 1 - win // 2, end=int(r["hg38_pos"]) + win // 2)
            if (wdir / "request.json").exists():
                wres = aa.load_archive(wdir)
            else:
                t0 = time.time()
                wres = aq.call_with_quota_retry(lambda: client.query_interval(interval, requested_scorers=SCORERS, progress_bar=False, max_workers=4), label=r["rsid"])
                aa.archive_scores(wres, wdir, {"rsid": r["rsid"], "interval": f"{r['hg38_chrom']}:{interval.start + 1}-{interval.end}", "requested_scorers": SCORERS,
                                               "sdk_version": sdk, "elapsed_seconds": time.time() - t0})
            if r["served"]:
                for res in point.values():
                    v, e = channel_summary(res, r["variant_uid"], r["ensembl"])
                    if any(x == x for x in v.values()):
                        vals, extra = v, e
                # background: every other SNV in the window, same channel rule
                bg = defaultdict(list)
                uids_in_window = set()
                for scorer_res in [wres]:
                    for ch, scorer in CHANNELS.items():
                        a = scorer_res.get(scorer)
                        if a is None or a.shape[0] == 0 or "variant" not in a.obs:
                            continue
                        obs = a.obs.reset_index(drop=True)
                        mask = liver_mask(a.var, "H3K27ac" if ch == "h3k27ac" else None)
                        if not mask.any():
                            continue
                        q = np.asarray(a.layers["quantiles"], np.float32)[:, mask]
                        for i in range(len(obs)):
                            u = aa.variant_uid_from_str(str(obs.at[i, "variant"]))
                            uids_in_window.add(u)
                            if u == r["variant_uid"]:
                                continue
                            if scorer in ("RNA_SEQ", "SPLICE_SITE_USAGE") and str(obs.at[i, "gene_id"]).split(".")[0] != r["ensembl"]:
                                continue
                            bg[ch].append(float(np.nanmax(np.abs(q[i]))) if ch == "splice_site_usage" else float(np.nanmedian(q[i])))
                ranks = {ch: rank_percentile(vals[ch], bg[ch]) for ch in CHANNELS}
                extra["n_window_snvs"] = len(uids_in_window)
            else:
                extra["window_covers_position"] = any(
                    int(aa.variant_uid_from_str(str(v)).split(":")[1]) == int(r["hg38_pos"])
                    for a in wres.values() if a.shape[0] and "variant" in a.obs for v in a.obs["variant"])
        returned = any(v == v for v in vals.values())
        verdict = score_expectation(r["expected_channel"], r["served"], returned, vals, ranks)
        out_rows.append({**{k: r[k] for k in ("rsid", "gene", "ensembl", "class", "expected_channel", "risk_allele", "variant_uid", "hg38_chrom", "hg38_pos", "ref", "alt", "allele_state", "served")}, "atlas_returned_a_row": returned,
                         **{f"{c}_liver_quantile": vals[c] for c in CHANNELS}, **{f"{c}_window_rank": ranks[c] for c in CHANNELS}, **extra, "verdict": verdict, "note": r["note"]})
    cols = list(dict.fromkeys(list(out_rows[0].keys()) + sorted({k for r in out_rows for k in r})))
    la.write_tsv_once(TABLES / "anchor_panel.tsv", out_rows, cols)
    coding = [r for r in out_rows if r["class"].startswith("coding_missense")]
    summary = {
        "n_anchors": len(out_rows), "n_served": sum(r["served"] for r in out_rows), "verdicts": {v: sum(r["verdict"] == v for r in out_rows) for v in ("hit", "miss", "uncovered", "not_applicable", "not_served", "not_prespecified")},
        "predictions_scored": {
            "P2.1": all(abs(r["rna_gene_liver_quantile"]) < 0.5 for r in coding if r["rna_gene_liver_quantile"] == r["rna_gene_liver_quantile"]),
            "P2.2": next((r["allele_state"] == "indel_not_served_by_point_endpoint" for r in out_rows if r["rsid"] == "rs72613567"), None),
            "P2.3": next(((r["rna_gene_liver_quantile"] == r["rna_gene_liver_quantile"] and abs(r["rna_gene_liver_quantile"]) >= 0.3 and
                           ((r["alt"] == "T") == (r["rna_gene_liver_quantile"] < 0))) or any((r[f"{c}_window_rank"] or 0) >= 0.9 for c in ("atac", "dnase", "h3k27ac"))
                          for r in out_rows if r["rsid"] == "rs641738"), None),
            "P2.4": sum(any((r[f"{c}_window_rank"] == r[f"{c}_window_rank"]) and r[f"{c}_window_rank"] >= 0.9 for c in ("atac", "dnase", "h3k27ac")) for r in coding) <= 2,
        },
        "per_anchor": {r["rsid"]: {"verdict": r["verdict"], "rna_gene": r["rna_gene_liver_quantile"], "avi_quantile": r.get("avi_quantile")} for r in out_rows},
    }
    json.dump(summary, (TABLES / "anchor_panel_summary.json").open("w"), indent=1, default=float)
    la.log(f"anchor panel: {summary['verdicts']}; predictions {summary['predictions_scored']}")


if __name__ == "__main__":
    main()
