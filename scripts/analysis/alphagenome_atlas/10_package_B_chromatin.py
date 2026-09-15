#!/usr/bin/env python3
"""Step 10 (Package B): predicted regulatory susceptibility beside measured liver chromatin.

Census only. Static overlap across five lineages; disease change restricted to hepatocyte and
stellate (the only lineages with jointly testable peaks). Posterior mass enters once per
variant per lineage (rule copied from atac_context_v3/09). Matched-background comparison is
descriptive: full control distribution, no headline p-value.

Outputs (tables/): chromatin_census.tsv, variant_peak_overlaps.tsv.gz, measured_context.tsv,
                   chromatin_matched_background.tsv
"""

from __future__ import annotations

import bisect
import csv
import math
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

import lib_atlas as la

P = la.prespec()
TABLES = la.out_root() / "tables"
CTX = la.PROJECT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
LINEAGES = ["hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"]
DA_LINEAGES = ["hepatocyte", "stellate"]
DRAWS, SEED = P["resampling"]["draws"], P["seeds"]["master"]
GENETIC_MIN_W = 0.01


class PeakIndex:
    """Sorted starts + prefix-max ends; returns the peak ids overlapping a 1-based position."""

    def __init__(self):
        self.items = defaultdict(list)
        self.starts, self.prefix, self.rows = {}, {}, {}

    def add(self, lineage, chrom, start0, end, pid):
        self.items[(lineage, chrom)].append((int(start0), int(end), pid))

    def build(self):
        for key, iv in self.items.items():
            iv.sort()
            self.rows[key] = iv
            self.starts[key] = [v[0] for v in iv]
            run, pre = -1, []
            for _, e, _ in iv:
                run = max(run, e); pre.append(run)
            self.prefix[key] = pre

    def hits(self, lineage, chrom, pos1):
        key = (lineage, chrom)
        if key not in self.starts:
            return []
        up = bisect.bisect_left(self.starts[key], pos1)
        if up == 0 or self.prefix[key][up - 1] < pos1:
            return []
        out = []
        for s0, e, pid in self.rows[key][:up][::-1]:
            if s0 < pos1 <= e:
                out.append(pid)
            if self.prefix[key][up - 1] < pos1:
                break
        return out


def main() -> None:
    signals = {s["signal_uid"]: s for s in la.read_tsv(TABLES / "eligible_signals.tsv")}
    # peaks
    manifest = la.read_tsv(CTX / "consensus_peak_manifest.tsv")
    idx = PeakIndex()
    peak_info = {}
    for r in manifest:
        pid = f"{r['lineage']}|{r['chrom']}:{r['start0']}-{r['end']}"
        peak_info[pid] = r
        idx.add(r["lineage"], r["chrom"], r["start0"], r["end"], pid)
    idx.build()
    da = pd.read_csv(CTX / "da/da_peak_results.tsv.gz", sep="\t", keep_default_na=False)
    da["pid"] = da["lineage"] + "|" + da["peak_coordinate"]
    da = da.set_index("pid")
    summary = la.read_tsv(CTX / "da/da_lineage_summary.tsv")
    la.write_tsv_once(TABLES / "measured_context.tsv", [
        {"cohort": c, "lineage": r["lineage"], "n_normal": r[f"n_normal_gse{g}"], "n_mash": r[f"n_mash_gse{g}"],
         "n_jointly_testable": r["n_jointly_testable"], "n_sig_this_cohort": r[f"n_sig_gse{g}"],
         "n_supported": r["n_supported"], "n_discordant": r["n_discordant"], "state": r["lineage_state"], "unit": "donor", "contrast": "MASH vs normal",
         "existing_family": "atac-context-v3 BH within jointly testable peaks per cohort"}
        for r in summary for c, g in (("GSE244832", "244832"), ("GSE281367", "281367"))],
        ["cohort", "lineage", "n_normal", "n_mash", "n_jointly_testable", "n_sig_this_cohort", "n_supported", "n_discordant", "state", "unit", "contrast", "existing_family"])

    # Atlas accessibility per variant (liver tracks)
    liver = pd.read_csv(TABLES / "liver_summaries.tsv.gz", sep="\t", keep_default_na=False)
    atac = liver[liver["scorer"] == "ATAC"].set_index("variant_uid")
    dnase = liver[liver["scorer"] == "DNASE"].set_index("variant_uid")

    # variants with weights (all mapped variants, not only queried)
    census = defaultdict(lambda: defaultdict(float))
    counts = defaultdict(lambda: defaultdict(int))
    overlap_rows = []
    genetic_peaks = defaultdict(lambda: defaultdict(float))  # lineage -> pid -> max weight
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["mapping_status"] != "mapped":
                continue
            sig, v, w = r["signal_uid"], r["variant_uid"], float(r["weight"])
            chrom, pos = v.split(":")[0], int(v.split(":")[1])
            census[sig]["mapped_mass"] += w
            queried = r["in_query_set"] == "True"
            a_raw = float(atac.at[v, "liver_median_raw"]) if v in atac.index else math.nan
            a_q = float(atac.at[v, "liver_median_quantile"]) if v in atac.index else math.nan
            d_raw = float(dnase.at[v, "liver_median_raw"]) if v in dnase.index else math.nan
            if not math.isnan(a_raw):
                census[sig]["atac_C"] += w; census[sig]["atac_A"] += w * abs(a_raw); census[sig]["atac_S"] += w * a_raw
                census[sig]["atac_A_quantile"] += w * abs(a_q)
            for lin in LINEAGES:
                hits = idx.hits(lin, chrom, pos)
                if not hits:
                    continue
                acc = {c: any(peak_info[p][f"supported_{c}"] == "TRUE" for p in hits) for c in ("gse244832", "gse281367")}
                both = any(peak_info[p]["supported_both"] == "TRUE" for p in hits)
                if acc["gse244832"] or acc["gse281367"]:
                    census[sig][f"{lin}_any_mass"] += w; counts[sig][f"{lin}_n_any"] += 1
                if both:
                    census[sig][f"{lin}_both_mass"] += w
                for c, ok in acc.items():
                    if ok:
                        census[sig][f"{lin}_{c}_mass"] += w
                if not math.isnan(a_raw) and (acc["gse244832"] or acc["gse281367"]):
                    census[sig][f"{lin}_inpeak_atac_C"] += w; census[sig][f"{lin}_inpeak_atac_A"] += w * abs(a_raw)
                if lin in DA_LINEAGES:
                    for p in hits:
                        if p not in da.index:
                            continue
                        d = da.loc[p]
                        state = d["evidence_state"]
                        census[sig][f"{lin}_da_{state}_mass"] += w
                        if w >= GENETIC_MIN_W:
                            genetic_peaks[lin][p] = max(genetic_peaks[lin][p], w)
                        overlap_rows.append({
                            "signal_uid": sig, "universe": signals[sig]["universe"], "variant_uid": v, "weight": w, "queried": queried, "lineage": lin,
                            "peak": p.split("|")[1], "peak_width": int(peak_info[p]["end"]) - int(peak_info[p]["start0"]),
                            "accessible_gse244832": acc["gse244832"], "accessible_gse281367": acc["gse281367"],
                            "logFC_gse244832": d["logFC_gse244832"], "se_gse244832": d["standard_error_gse244832"], "q_gse244832": d["qvalue_gse244832"],
                            "logFC_gse281367": d["logFC_gse281367"], "se_gse281367": d["standard_error_gse281367"], "q_gse281367": d["qvalue_gse281367"],
                            "evidence_state": state, "promoter_genes": d["promoter_genes"],
                            "atlas_atac_liver_median_raw": a_raw, "atlas_atac_liver_median_quantile": a_q, "atlas_dnase_liver_median_raw": d_raw,
                            "window_compatibility": "variant_inside_500bp_peak; Atlas ATAC scorer window not exposed by the API (center-mask, width unspecified)",
                            "sign_caveat": "allele effect sign is not expected to equal the MASH-vs-normal sign",
                        })
                        break  # one DA row per variant per lineage (peaks do not overlap within a lineage consensus set)
    census_rows = []
    for sig, s in signals.items():
        c = census[sig]
        row = {"signal_uid": sig, "universe": s["universe"], "gwas_name": s["gwas_name"], "trait_class": s["trait_class"], "gene": s["gene"],
               "mapped_mass": c["mapped_mass"], "atac_C": c["atac_C"], "atac_A_over_C": (c["atac_A"] / c["atac_C"]) if c["atac_C"] else math.nan,
               "atac_S_over_C": (c["atac_S"] / c["atac_C"]) if c["atac_C"] else math.nan,
               "atac_A_quantile_over_C": (c["atac_A_quantile"] / c["atac_C"]) if c["atac_C"] else math.nan}
        for lin in LINEAGES:
            row[f"{lin}_any_mass"] = c[f"{lin}_any_mass"]; row[f"{lin}_both_mass"] = c[f"{lin}_both_mass"]
            row[f"{lin}_gse244832_mass"] = c[f"{lin}_gse244832_mass"]; row[f"{lin}_gse281367_mass"] = c[f"{lin}_gse281367_mass"]
            row[f"{lin}_inpeak_atac_A_over_C"] = (c[f"{lin}_inpeak_atac_A"] / c[f"{lin}_inpeak_atac_C"]) if c[f"{lin}_inpeak_atac_C"] else math.nan
            if lin in DA_LINEAGES:
                for st in ("supported", "discordant", "source_dependent", "indeterminate"):
                    row[f"{lin}_da_{st}_mass"] = c[f"{lin}_da_{st}_mass"]
            else:
                row[f"{lin}_disease_change"] = "unavailable_no_jointly_testable_peaks" if lin == "macrophage" else "unavailable_not_in_da_design"
        census_rows.append(row)
    la.write_tsv_once(TABLES / "chromatin_census.tsv", census_rows, list(census_rows[0].keys()))
    la.write_tsv_once(TABLES / "variant_peak_overlaps.tsv.gz", overlap_rows, list(overlap_rows[0].keys()) if overlap_rows else ["signal_uid"])

    # matched background (descriptive): genetic peaks (w ≥ 0.01) vs draws of testable peaks matched on promoter status and GC decile
    fasta = pysam.FastaFile(FASTA)
    rng = np.random.default_rng(SEED)
    bg_rows = []
    for lin in DA_LINEAGES:
        d = da[da["lineage"] == lin].copy()
        d["gc"] = [gc_content(fasta, p) for p in d["peak_coordinate"]]
        d["promoter"] = d["promoter_genes"].astype(str) != ""
        d["gc_decile"] = pd.qcut(d["gc"], 10, labels=False, duplicates="drop")
        d["stratum"] = d["promoter"].astype(str) + "|" + d["gc_decile"].astype(str)
        gen = d[d.index.isin(genetic_peaks[lin])]
        if len(gen) == 0:
            bg_rows.append({"lineage": lin, "n_genetic_peaks": 0, "note": "no genetic peak at weight >= 0.01"})
            continue
        strata_counts = gen["stratum"].value_counts()
        pools = {s: d[(d["stratum"] == s) & (~d.index.isin(gen.index))] for s in strata_counts.index}
        for cohort in ("gse244832", "gse281367"):
            obs_abs = gen[f"logFC_{cohort}"].astype(float).abs().mean()
            obs_pos = (gen[f"logFC_{cohort}"].astype(float) > 0).mean()
            draws_abs, draws_pos = np.empty(DRAWS), np.empty(DRAWS)
            for k in range(DRAWS):
                parts = [pools[s].sample(n=min(n, len(pools[s])), random_state=int(rng.integers(0, 2**31 - 1))) for s, n in strata_counts.items() if len(pools[s])]
                m = pd.concat(parts)
                vals = m[f"logFC_{cohort}"].astype(float)
                draws_abs[k], draws_pos[k] = vals.abs().mean(), (vals > 0).mean()
            bg_rows.append({"lineage": lin, "cohort": cohort, "n_genetic_peaks": len(gen), "n_testable_peaks": len(d),
                            "matching": "promoter_status+GC_decile (fixed width 500 bp); baseline abundance/detection/mappability not available in the DA table",
                            "observed_mean_abs_logFC": obs_abs, "control_mean_abs_logFC_median": float(np.median(draws_abs)),
                            "control_abs_q025": float(np.quantile(draws_abs, 0.025)), "control_abs_q975": float(np.quantile(draws_abs, 0.975)),
                            "exceedances_abs": int((draws_abs >= obs_abs).sum()), "observed_fraction_positive": obs_pos,
                            "control_fraction_positive_median": float(np.median(draws_pos)), "draws": DRAWS, "seed": SEED,
                            "inference": "descriptive; matched draws do not preserve genomic dependence, no calibrated p-value claimed"})
    la.write_tsv_once(TABLES / "chromatin_matched_background.tsv", bg_rows, sorted({k for r in bg_rows for k in r}))
    la.log(f"Package B: {len(census_rows)} signals, {len(overlap_rows)} variant-peak rows, background rows {len(bg_rows)}")


def gc_content(fasta, coord: str) -> float:
    chrom, rng_ = coord.split(":")
    s0, e = map(int, rng_.split("-"))
    seq = fasta.fetch(chrom, s0, e).upper()
    return (seq.count("G") + seq.count("C")) / max(len(seq), 1)


if __name__ == "__main__":
    main()
