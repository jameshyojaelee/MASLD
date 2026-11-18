#!/usr/bin/env python3
"""DIG1 62b (2026-11-18): Genotype calling with strict HWE >= 0.05, dp>=20 sensitivity, and autosomal ancestry PCs."""
from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
DIG1 = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
RNA_GT = DIG1 / "out_rna_genotypes"
ANC_TSV = ROOT / "RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else DIG1 / "out_genotype_table_20261118"
OUT.mkdir(parents=True, exist_ok=True)


def hwe_exact(obs_hom1: int, obs_hets: int, obs_hom2: int) -> float:
    """Wigginton exact two-sided Hardy-Weinberg test p-value."""
    if obs_hom1 < 0 or obs_hets < 0 or obs_hom2 < 0:
        return float("nan")
    n = obs_hom1 + obs_hets + obs_hom2
    if n == 0:
        return float("nan")
    obs_homr = min(obs_hom1, obs_hom2)
    rare_copies = 2 * obs_homr + obs_hets
    if rare_copies == 0:
        return 1.0
    probs = [0.0] * (rare_copies + 1)
    mid = (rare_copies * (2 * n - rare_copies)) // (2 * n)
    if (rare_copies & 1) != (mid & 1):
        mid += 1
    curr_hets = mid
    curr_homr = (rare_copies - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    probs[mid] = 1.0
    tot = 1.0
    while curr_hets > 1:
        p = probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        probs[curr_hets - 2] = p
        tot += p
        curr_hets -= 2
        curr_homr += 1
        curr_homc += 1
    curr_hets = mid
    curr_homr = (rare_copies - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets <= rare_copies - 2:
        p = probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        probs[curr_hets + 2] = p
        tot += p
        curr_hets += 2
        curr_homr -= 1
        curr_homc -= 1
    target = probs[obs_hets]
    return min(1.0, sum(p for p in probs if p <= target + 1e-12) / tot)


def gt_to_dosage(gt_str: str) -> float:
    gt = gt_str.split(":", 1)[0]
    if gt in {"0/0", "0|0"}:
        return 0.0
    if gt in {"0/1", "1/0", "0|1", "1|0"}:
        return 1.0
    if gt in {"1/1", "1|1"}:
        return 2.0
    return float("nan")


def call_dosage(depth: int, alt_frac: float, min_dp: int) -> float:
    if depth < min_dp or not np.isfinite(alt_frac):
        return float("nan")
    if alt_frac <= 0.10:
        return 0.0
    if 0.25 <= alt_frac <= 0.75:
        return 1.0
    if alt_frac >= 0.90:
        return 2.0
    return float("nan")


def parse_vcf_sample(vcf_path: Path, target_map: dict[tuple[str, int], tuple[str, str]]) -> list[dict]:
    records = []
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _vid, _ref, _alt, _qual, _flt, info, fmt, sample_val = parts[:10]
            key = (chrom, int(pos))
            if key not in target_map:
                continue
            rsid, gene = target_map[key]
            fmt_keys = fmt.split(":")
            fmt_vals = sample_val.split(":")
            fdict = dict(zip(fmt_keys, fmt_vals))
            gt_bcf = gt_to_dosage(fdict.get("GT", "./."))
            ad_str = fdict.get("AD", "")
            ref_n, alt_n = 0, 0
            if ad_str and ad_str != ".":
                ad_parts = [int(x) for x in ad_str.split(",") if x != "."]
                if len(ad_parts) >= 1:
                    ref_n = ad_parts[0]
                if len(ad_parts) >= 2:
                    alt_n = sum(ad_parts[1:])
            else:
                # Fallback to INFO/DP4
                for kv in info.split(";"):
                    if kv.startswith("DP4="):
                        d4 = [int(x) for x in kv[4:].split(",")]
                        ref_n = d4[0] + d4[1]
                        alt_n = d4[2] + d4[3]
            depth = ref_n + alt_n
            alt_frac = (alt_n / float(depth)) if depth > 0 else float("nan")
            records.append({
                "rsid": rsid,
                "gene": gene,
                "depth": int(depth),
                "ref_n": int(ref_n),
                "alt_n": int(alt_n),
                "alt_frac": alt_frac,
                "call": call_dosage(depth, alt_frac, 10),
                "call_dp20": call_dosage(depth, alt_frac, 20),
                "gt_bcftools": gt_bcf,
            })
    return records


def cohort_qc_stats(sub: pd.DataFrame, call_col: str) -> dict:
    n_tot = len(sub)
    called = sub[sub[call_col].notna()]
    n_called = len(called)
    cr = n_called / float(n_tot) if n_tot > 0 else float("nan")
    n0 = int((called[call_col] == 0.0).sum())
    n1 = int((called[call_col] == 1.0).sum())
    n2 = int((called[call_col] == 2.0).sum())
    af = (n1 + 2.0 * n2) / (2.0 * n_called) if n_called > 0 else float("nan")
    hw = hwe_exact(n0, n1, n2) if n_called > 0 else float("nan")
    both = called[called["gt_bcftools"].notna()]
    agree = float((both[call_col] == both["gt_bcftools"]).mean()) if len(both) > 0 else float("nan")
    return dict(n=n_tot, n_called=n_called, call_rate=cr, n0=n0, n1=n1, n2=n2,
                alt_freq=af, hwe_p=hw, agree_bcftools=agree)


def main():
    samples = pd.read_csv(RNA_GT / "samples.tsv", sep="\t")
    targets = pd.read_csv(
        RNA_GT / "62_targets_annotated.tsv", sep="\t", header=None,
        names=["chrom", "start", "end", "rsid", "gene", "ref", "alt", "note"],
    )
    target_map = {(r.chrom, int(r.start)): (r.rsid, r.gene) for r in targets.itertuples(index=False)}
    all_rsids = targets["rsid"].tolist()

    anc = pd.read_csv(ANC_TSV, sep="\t")
    anc_cols = ["run", "nearest_superpop", "ancestry_call", "passes_min_sites", "cohort_contaminated"] + [f"PC{i}" for i in range(1, 11)]
    anc_sub = anc[anc_cols].rename(columns={"run": "sample_id"})
    samples = samples.merge(anc_sub, on="sample_id", how="left")
    assert samples["PC1"].notna().all(), "Missing ancestry PCs for some samples!"

    rows = []
    vcf_dir = RNA_GT / "vcf"
    for r in samples.itertuples(index=False):
        vp = vcf_dir / f"{r.sample_id}.vcf.gz"
        recs = parse_vcf_sample(vp, target_map)
        rec_by_rs = {x["rsid"]: x for x in recs}
        for t in targets.itertuples(index=False):
            x = rec_by_rs.get(t.rsid, {
                "rsid": t.rsid, "gene": t.gene, "depth": 0, "ref_n": 0, "alt_n": 0,
                "alt_frac": float("nan"), "call": float("nan"), "call_dp20": float("nan"),
                "gt_bcftools": float("nan"),
            })
            rows.append({
                "sample_id": r.sample_id,
                "dataset": r.dataset,
                "nearest_superpop": r.nearest_superpop,
                "ancestry_call": r.ancestry_call,
                **x,
            })

    long_df = pd.DataFrame(rows)
    long_df.to_csv(OUT / "genotype_long.tsv.gz", sep="\t", index=False, compression="gzip")

    # Verify exact concordance with prior dp>=10 genotype_long.tsv.gz
    orig_path = DIG1 / "out_genotype_table/genotype_long.tsv.gz"
    if orig_path.exists():
        orig = pd.read_csv(orig_path, sep="\t")
        m_chk = long_df.merge(orig[["sample_id", "rsid", "call", "depth"]], on=["sample_id", "rsid"], suffixes=("_new", "_old"))
        assert len(m_chk) == len(orig), "Row count mismatch vs original genotype_long!"
        same_na = (m_chk["call_new"].isna() == m_chk["call_old"].isna()).all()
        same_val = (m_chk.loc[m_chk["call_new"].notna(), "call_new"] == m_chk.loc[m_chk["call_new"].notna(), "call_old"]).all()
        assert same_na and same_val, "Parsed call values differ from original 62b genotype_long!"

    wide_dp10 = long_df.pivot(index=["sample_id", "dataset"], columns="rsid", values="call")[all_rsids].reset_index()
    wide_dp20 = long_df.pivot(index=["sample_id", "dataset"], columns="rsid", values="call_dp20")[all_rsids].reset_index()
    wide_dp10.to_csv(OUT / "genotype_calls_wide.tsv", sep="\t", index=False)
    wide_dp20.to_csv(OUT / "genotype_calls_wide_dp20.tsv", sep="\t", index=False)

    # Cohort-by-site QC across dp>=10, dp>=20, and EUR-only subsets
    by_cohort_rows = []
    for (rsid, ds), sub in long_df.groupby(["rsid", "dataset"]):
        s10 = cohort_qc_stats(sub, "call")
        s20 = cohort_qc_stats(sub, "call_dp20")
        sub_eur = sub[sub["nearest_superpop"] == "EUR"]
        se10 = cohort_qc_stats(sub_eur, "call")
        se20 = cohort_qc_stats(sub_eur, "call_dp20")
        n_amb10 = int(((sub["depth"] >= 10) & sub["call"].isna()).sum())
        n_amb20 = int(((sub["depth"] >= 20) & sub["call_dp20"].isna()).sum())
        by_cohort_rows.append({
            "rsid": rsid,
            "cohort": ds,
            "n": s10["n"],
            "n_called": s10["n_called"],
            "call_rate": s10["call_rate"],
            "median_depth": float(sub["depth"].median()),
            "n0": s10["n0"],
            "n1": s10["n1"],
            "n2": s10["n2"],
            "alt_freq": s10["alt_freq"],
            "hwe_p": s10["hwe_p"],
            "agree_bcftools": s10["agree_bcftools"],
            "n_ambiguous": n_amb10,
            "n_called_dp20": s20["n_called"],
            "call_rate_dp20": s20["call_rate"],
            "n0_dp20": s20["n0"],
            "n1_dp20": s20["n1"],
            "n2_dp20": s20["n2"],
            "alt_freq_dp20": s20["alt_freq"],
            "hwe_p_dp20": s20["hwe_p"],
            "n_ambiguous_dp20": n_amb20,
            "n_eur_called": se10["n_called"],
            "n0_eur": se10["n0"],
            "n1_eur": se10["n1"],
            "n2_eur": se10["n2"],
            "alt_freq_eur": se10["alt_freq"],
            "hwe_p_eur": se10["hwe_p"],
            "n_eur_called_dp20": se20["n_called"],
            "hwe_p_eur_dp20": se20["hwe_p"],
        })
    per_cohort = pd.DataFrame(by_cohort_rows).sort_values(["rsid", "cohort"])
    per_cohort.to_csv(OUT / "site_qc_by_cohort.tsv", sep="\t", index=False)

    # Aggregate site QC with STRICT per-cohort hwe_p >= 0.05 (among cohorts with >= 20 calls)
    site_rows = []
    for rsid in all_rsids:
        sub = long_df[long_df["rsid"] == rsid]
        per = per_cohort[per_cohort["rsid"] == rsid]
        s10 = cohort_qc_stats(sub, "call")
        s20 = cohort_qc_stats(sub, "call_dp20")
        sub_eur = sub[sub["nearest_superpop"] == "EUR"]
        se10 = cohort_qc_stats(sub_eur, "call")

        hw10_vals = per.loc[per["n_called"] >= 20, "hwe_p"].dropna()
        min_hw10 = float(hw10_vals.min()) if len(hw10_vals) else float("nan")
        hw20_vals = per.loc[per["n_called_dp20"] >= 20, "hwe_p_dp20"].dropna()
        min_hw20 = float(hw20_vals.min()) if len(hw20_vals) else float("nan")
        hweur10_vals = per.loc[per["n_eur_called"] >= 20, "hwe_p_eur"].dropna()
        min_hweur10 = float(hweur10_vals.min()) if len(hweur10_vals) else float("nan")
        hweur20_vals = per.loc[per["n_eur_called_dp20"] >= 20, "hwe_p_eur_dp20"].dropna()
        min_hweur20 = float(hweur20_vals.min()) if len(hweur20_vals) else float("nan")

        min_cr10 = float(per["call_rate"].min())
        min_cr20 = float(per["call_rate_dp20"].min())
        se10_cohort_cr = (per["n_eur_called"] / np.maximum(1, long_df[long_df["rsid"] == rsid].groupby("dataset")["nearest_superpop"].apply(lambda x: (x == "EUR").sum()).reindex(per["cohort"]).to_numpy())).min()
        min_creur20 = (per["n_eur_called_dp20"] / np.maximum(1, long_df[long_df["rsid"] == rsid].groupby("dataset")["nearest_superpop"].apply(lambda x: (x == "EUR").sum()).reindex(per["cohort"]).to_numpy())).min()

        pass_strict = bool(
            (s10["call_rate"] >= 0.80)
            and (min_cr10 >= 0.50)
            and np.isfinite(s10["agree_bcftools"]) and (s10["agree_bcftools"] >= 0.98)
            and np.isfinite(min_hw10) and (min_hw10 >= 0.05)
        )
        pass_dp20 = bool(
            (s20["call_rate"] >= 0.80)
            and (min_cr20 >= 0.50)
            and np.isfinite(s20["agree_bcftools"]) and (s20["agree_bcftools"] >= 0.98)
            and np.isfinite(min_hw20) and (min_hw20 >= 0.05)
        )
        pass_eur05 = bool(
            (se10["call_rate"] >= 0.80)
            and (float(se10_cohort_cr) >= 0.50)
            and np.isfinite(se10["agree_bcftools"]) and (se10["agree_bcftools"] >= 0.98)
            and np.isfinite(min_hweur10) and (min_hweur10 >= 0.05)
        )
        pass_eur20_05 = bool(
            (cohort_qc_stats(sub_eur, "call_dp20")["call_rate"] >= 0.80)
            and (float(min_creur20) >= 0.50)
            and np.isfinite(min_hweur20) and (min_hweur20 >= 0.05)
        )
        pass_legacy_1e4 = bool(
            (s10["call_rate"] >= 0.80)
            and (min_cr10 >= 0.50)
            and np.isfinite(s10["agree_bcftools"]) and (s10["agree_bcftools"] >= 0.98)
            and np.isfinite(min_hw10) and (min_hw10 >= 1e-4)
        )
        site_rows.append({
            "rsid": rsid,
            "call_rate": s10["call_rate"],
            "min_cohort_call_rate": min_cr10,
            "min_cohort_hwe_p": min_hw10,
            "agree_bcftools": s10["agree_bcftools"],
            "alt_freq": s10["alt_freq"],
            "pass_qc": pass_strict,
            "call_rate_dp20": s20["call_rate"],
            "min_cohort_hwe_p_dp20": min_hw20,
            "agree_bcftools_dp20": s20["agree_bcftools"],
            "alt_freq_dp20": s20["alt_freq"],
            "pass_qc_dp20": pass_dp20,
            "call_rate_eur": se10["call_rate"],
            "min_cohort_hwe_p_eur": min_hweur10,
            "min_cohort_hwe_p_eur_dp20": min_hweur20,
            "pass_qc_eur_hwe05": pass_eur05,
            "pass_qc_eur_dp20_hwe05": pass_eur20_05,
            "pass_qc_hwe_1e4_legacy": pass_legacy_1e4,
        })
    site_qc = pd.DataFrame(site_rows)
    site_qc.to_csv(OUT / "site_qc.tsv", sep="\t", index=False)

    pip_out = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
    (OUT / "pip_freeze.txt").write_text(pip_out)
    print(site_qc[["rsid", "call_rate", "min_cohort_hwe_p", "pass_qc", "min_cohort_hwe_p_dp20", "pass_qc_dp20", "min_cohort_hwe_p_eur", "pass_qc_eur_hwe05", "pass_qc_hwe_1e4_legacy"]].to_string(index=False))


if __name__ == "__main__":
    main()
