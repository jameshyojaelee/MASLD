#!/usr/bin/env python3
"""Step 33 (P3a): evaluate Atlas scores on the external same-variant benchmarks (apply-only replications).

usage: 33_benchmark_eval.py caqtl|sqtl|mpra|mpra_rescore

mpra_rescore re-scores from the deposited per-interval table without re-reading the interval archives
(176 GB, about 40 minutes), which is what you want after a scoring-side fix.

caqtl : Currin 2025 bulk-liver caQTL leads. The estimator (auROC, chromosome-block bootstrap seed 123, label
        permutation seed 99, 500-kb signal clump, positional baseline excluded here) is copied from
        GWAS/finemapping/src/84_caqtl_zeroshot.py so the Atlas ATAC row can be compared with the archived
        AlphaGenome-API row (0.707 variant level / 0.726 signal level). One row per Atlas channel; nothing tuned.
sqtl  : GTEx v8 liver sGene leads and per-cluster top pairs (positives only): sign concordance of slope with the
        Atlas splice-site-usage / junction liver quantile at the matching intron, and Spearman(|slope|, |quantile|).
mpra  : Hu/Zhu 2026 MPRA oligos: DAV vs non-DAV per context using the allele-agnostic centre-position maximum
        |quantile| per channel (rule fixed in 32_mpra_library.py), plus allelic sign where alleles are known.

Outputs (tables/): external_benchmark_<set>.tsv, external_benchmark_<set>.json
"""

from __future__ import annotations

import json
import math
import sys
import collections
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

import atlas_archive as aa
import lib_atlas as la

ROOT = la.out_root()
RAW = ROOT / "raw"
TABLES = ROOT / "tables"
CAQTL = la.PROJECT / "GWAS/finemapping/results/seqfunc/direction_features_caqtl.tsv"
PRIOR = la.PROJECT / "GWAS/finemapping/results/seqfunc/caqtl_zeroshot_permodel.tsv"
SQTL_DIR = la.PROJECT / "data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL"
LIVER = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
CLUMP_BP = 500_000
# Fixed before retrieval and not adjusted after any MPRA score was seen (the AlphaGenome exposure audit
# forbids MPRA outcomes selecting tracks or aggregation). One definition, so the read path and the re-score
# path cannot drift apart. fibroblast_dnase is a CROSS-TISSUE proxy for LX-2; the Atlas has no stellate track.
MPRA_CHANNELS = {"hepg2_atac": ("ATAC", "HepG2"), "hepg2_h3k27ac": ("CHIP_HISTONE", "HepG2.*H3K27ac|H3K27ac.*HepG2"),
                 "hepg2_cage": ("CAGE", "HepG2"), "hepg2_dnase": ("DNASE", "HepG2"),
                 "liver_atac": ("ATAC", None), "liver_h3k27ac": ("CHIP_HISTONE", "H3K27ac"),
                 "fibroblast_dnase": ("DNASE", "fibroblast"), "avi": ("AVI_SCORE", None)}
NBOOT = 2000
NPERM = 1000
CAQTL_CHANNELS = {  # scorer -> (track filter on name, signed?)
    "atlas_atac": ("ATAC", None), "atlas_dnase": ("DNASE", None), "atlas_h3k27ac": ("CHIP_HISTONE", "H3K27ac"), "atlas_h3k4me1": ("CHIP_HISTONE", "H3K4me1"),
    "atlas_h3k4me3": ("CHIP_HISTONE", "H3K4me3"), "atlas_tf_hnf4a": ("CHIP_TF", "HNF4A"), "atlas_tf_foxa2": ("CHIP_TF", "FOXA2"), "atlas_tf_cebpb_hepg2": ("CHIP_TF", "CEBPB"),
    "atlas_avi": ("AVI_SCORE", None),
}


# ---------------------------------------------------------------- estimator, copied from 84_caqtl_zeroshot.py (sklearn-free auROC)
def auroc(y: np.ndarray, s: np.ndarray, ok=None) -> float:
    if ok is None:
        ok = ~np.isnan(s)
    if ok.sum() == 0 or len(np.unique(y[ok])) < 2:
        return np.nan
    yy, ss = y[ok], s[ok]
    r = rankdata(ss)                                   # average ranks = roc_auc_score tie handling
    n1 = int(yy.sum()); n0 = len(yy) - n1
    return float((r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def block_bootstrap(y, score, groups, ok, nboot=NBOOT, seed=123):
    rng = np.random.default_rng(seed)
    yo, so, go = y[ok], score[ok], np.asarray(groups)[ok]
    chroms = pd.unique(go)
    # A context/channel can select zero eligible rows, or rows of a single class. Neither supports an
    # interval; returning NaN keeps one empty cell from killing the whole benchmark, which is what happened
    # after the MPRA run had already read 176 GB of archives.
    if len(chroms) == 0 or len(np.unique(yo)) < 2:
        return (np.nan, np.nan)
    aucs = []
    for _ in range(nboot):
        pick = rng.choice(chroms, size=len(chroms), replace=True)
        idx = np.concatenate([np.where(go == c)[0] for c in pick])
        if len(np.unique(yo[idx])) < 2:
            continue
        aucs.append(auroc(yo[idx], so[idx], np.ones(len(idx), bool)))
    if not aucs:
        return (np.nan, np.nan)
    return (float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)))


def perm_null_p(y, score, ok, obs, nperm=NPERM, seed=99):
    rng = np.random.default_rng(seed)
    yo, so = y[ok], score[ok]
    null = []
    for _ in range(nperm):
        yp = rng.permutation(yo)
        if len(np.unique(yp)) < 2:
            continue
        null.append(auroc(yp, so, np.ones(len(yp), bool)))
    null = np.asarray(null)
    a = max(obs, 1 - obs)
    hit = np.sum((null >= a) | (null <= 1 - a))
    return float((1 + int(hit)) / (1 + len(null))), float(np.mean(null)) if len(null) else np.nan, len(null)


def signal_clump(df: pd.DataFrame) -> pd.DataFrame:
    keep = []
    for ch, sub in df.groupby("chr"):
        sub = sub.sort_values("pos_hg38").reset_index()
        last = -10 * CLUMP_BP
        cur, clumps = [], []
        for _, r in sub.iterrows():
            if r["pos_hg38"] - last > CLUMP_BP and cur:
                clumps.append(cur); cur = []
            cur.append(r); last = r["pos_hg38"]
        if cur:
            clumps.append(cur)
        for cl in clumps:
            cldf = pd.DataFrame(cl)
            rep = cldf.loc[pd.to_numeric(cldf["pvalue"], errors="coerce").idxmin()]
            keep.append(int(rep["index"]))
    return df.loc[keep].copy()


def scorer_is_signed() -> dict[str, bool]:
    """Server-declared signedness per scorer, read from the Track 0 deposit.

    Added 2026-09-15 after the sQTL arm was found to have run a sign concordance on SPLICE_JUNCTIONS, which
    the server declares unsigned and whose every quantile is positive. Any statistic that reads a SIGN must
    consult this; a magnitude statistic need not.
    """
    path = la.track0_root() / "tables" / "scorer_metadata.tsv"
    return {r["scorer"]: str(r["is_signed"]) == "True" for r in la.read_tsv(path)}


def sign_note(scorer: str, signed: dict[str, bool]) -> str:
    """Empty for a signed scorer; otherwise the reason this row is not a direction test."""
    if signed.get(scorer, False):
        return ""
    return (f"NOT A DIRECTION TEST: {scorer} is server-declared unsigned (is_signed=False), so the sign of "
            "its value is not a reference-to-alternate direction. Read the magnitude rows instead.")


def score_one(y, score, groups, ok, label) -> dict:
    a = auroc(y, score, ok)
    lo, hi = block_bootstrap(y, score, groups, ok)
    p, null_mean, n_null = perm_null_p(y, score, ok, a)
    return dict(model=label, auroc=a, ci_lo=lo, ci_hi=hi, perm_p=p, perm_null_mean=null_mean, n_perm=n_null, n=int(ok.sum()))


# ---------------------------------------------------------------- archive readers
def liver_track_mask(var: pd.DataFrame, name_filter: str | None) -> np.ndarray:
    m = var["ontology_curie"].isin(LIVER).values if "ontology_curie" in var else np.zeros(len(var), bool)
    if name_filter is not None and "name" in var:
        m = var["name"].astype(str).str.contains(name_filter).values & (m | var["biosample_name"].astype(str).str.contains("HepG2", case=False).values)
    return m


def collect_point_scores(stage_dir, channels: dict) -> dict[str, dict[str, float]]:
    """{channel: {variant_uid: liver-median signed quantile}} (AVI: the raw score)."""
    out = {c: {} for c in channels}
    for chunk in sorted(stage_dir.glob("chunk_*")):
        res = aa.load_archive(chunk)
        for label, (scorer, name_filter) in channels.items():
            a = res.get(scorer)
            if a is None or a.shape[0] == 0 or "variant" not in a.obs:
                continue
            obs = a.obs.reset_index(drop=True)
            if scorer == "AVI_SCORE":
                x = np.asarray(a.X, np.float32)[:, 0]
                for i in range(len(obs)):
                    out[label][aa.variant_uid_from_str(str(obs.at[i, "variant"]))] = float(x[i])
                continue
            mask = liver_track_mask(a.var.reset_index(drop=True), name_filter)
            if not mask.any():
                continue
            q = np.asarray(a.layers["quantiles"], np.float32)[:, mask]
            for i in range(len(obs)):
                out[label][aa.variant_uid_from_str(str(obs.at[i, "variant"]))] = float(np.nanmedian(q[i]))
    return out


# ---------------------------------------------------------------- benchmarks
def run_caqtl() -> None:
    signed = scorer_is_signed()
    df = pd.read_csv(CAQTL, sep="\t", dtype=str)
    df["chr"] = df["chr"].astype(str).str.replace("^chr", "", regex=True)
    df["pos_hg38"] = pd.to_numeric(df["pos_hg38"], errors="coerce").astype(int)
    df["variant_uid"] = "chr" + df["chr"] + ":" + df["pos_hg38"].astype(str) + ":" + df["ref"] + ":" + df["alt"]
    df["orient"] = pd.to_numeric(df["ag_orient_mult"], errors="coerce").fillna(1.0)
    scores = collect_point_scores(RAW / "atlas_caqtl", CAQTL_CHANNELS)
    for label in CAQTL_CHANNELS:
        df[label] = df["variant_uid"].map(scores[label]) * df["orient"]          # label frame, as 79/84 did
    df["ag_api_atac"] = pd.to_numeric(df["ag_atac_delta_labelframe"], errors="coerce")
    df["ag_api_dnase"] = pd.to_numeric(df["ag_dnase_delta_labelframe"], errors="coerce")
    prior = {(r["level"], r["model"]): r for r in la.read_tsv(PRIOR)}
    rows = []
    for level, d in (("variant_level_all_peakleads", df), ("signal_level_clumped", signal_clump(df))):
        y = (pd.to_numeric(d["sign"], errors="coerce").to_numpy() > 0).astype(int)
        groups = d["chr"].astype(str).to_numpy()
        for label in list(CAQTL_CHANNELS) + ["ag_api_atac", "ag_api_dnase"]:
            s = pd.to_numeric(d[label], errors="coerce").to_numpy()
            ok = ~np.isnan(s)
            r = score_one(y, s, groups, ok, label); r["level"] = level; r["n_pos"] = int(y[ok].sum())
            # the label here is a SIGN, so an unsigned channel's auROC is not a direction test
            r["sign_validity"] = ("" if label.startswith("ag_api_")     # the archived API deltas are signed
                                  else sign_note(CAQTL_CHANNELS[label][0], signed))
            key = (level, {"ag_api_atac": "alphagenome_atac", "ag_api_dnase": "alphagenome_dnase"}.get(label, ""))
            if key in prior:
                r["archived_api_auroc"] = float(prior[key]["auroc"]); r["archived_api_n"] = int(prior[key]["n"])
            rows.append(r)
        # magnitude among strong caQTLs: |beta| vs |quantile|
        eff = pd.to_numeric(d["effect"], errors="coerce").abs().to_numpy()
        for label in CAQTL_CHANNELS:
            s = pd.to_numeric(d[label], errors="coerce").abs().to_numpy(); ok = ~np.isnan(s) & ~np.isnan(eff)
            rows.append({"level": level, "model": f"{label}|magnitude", "spearman_abs_effect_vs_abs_score": float(spearmanr(eff[ok], s[ok]).correlation) if ok.sum() > 10 else math.nan, "n": int(ok.sum())})
    la.write_tsv_once(TABLES / "external_benchmark_caqtl.tsv", rows, sorted({k for r in rows for k in r}))
    api_full = [r for r in rows if r["model"] == "ag_api_atac" and r["level"] == "variant_level_all_peakleads"][0]
    atlas_full = [r for r in rows if r["model"] == "atlas_atac" and r["level"] == "variant_level_all_peakleads"][0]
    shared = df.dropna(subset=["atlas_atac", "ag_api_atac"])
    corr = float(spearmanr(shared["atlas_atac"], shared["ag_api_atac"]).correlation) if len(shared) > 10 else math.nan
    verdict = {
        "prior_archived": {"alphagenome_atac_variant_level": prior.get(("variant_level_all_peakleads", "alphagenome_atac"), {}).get("auroc"),
                            "alphagenome_atac_signal_level": prior.get(("signal_level_clumped", "alphagenome_atac"), {}).get("auroc")},
        "atlas_atac": {k: atlas_full[k] for k in ("auroc", "ci_lo", "ci_hi", "perm_p", "n")},
        "api_atac_recomputed_here_all_leads": {k: api_full[k] for k in ("auroc", "ci_lo", "ci_hi", "n")},
        "spearman_atlas_vs_api_atac_scores": corr, "n_shared": int(len(shared)),
        "reproduction_guard": "the archived API number used a 3,000-lead subsample (79_score_panels.py subsample rule); the Atlas is scored on every served lead and on the same subsample below",
        "note": "quantile-layer liver-track median (Atlas) vs raw delta mean (API): same model, different summary; sign agreement is the replication target",
    }
    # exact same subsample as 79: is_peak_lead first, then p asc, head 3000
    d = df.copy(); d["_p"] = pd.to_numeric(d["pvalue"], errors="coerce").fillna(1.0)
    d = d.sort_values(["_p"], kind="mergesort").head(3000)
    y = (pd.to_numeric(d["sign"], errors="coerce").to_numpy() > 0).astype(int); groups = d["chr"].astype(str).to_numpy()
    for label in ("atlas_atac", "ag_api_atac"):
        s = pd.to_numeric(d[label], errors="coerce").to_numpy(); ok = ~np.isnan(s)
        verdict[f"{label}_on_79_subsample_3000"] = score_one(y, s, groups, ok, label)
    arch = verdict["prior_archived"]["alphagenome_atac_variant_level"]
    arch_n = prior.get(("variant_level_all_peakleads", "alphagenome_atac"), {}).get("n")
    verdict["G_api_reproduces_archive"] = reproduction_guard(
        recomputed_auroc=api_full["auroc"], recomputed_n=api_full["n"],
        archived_auroc=float(arch), archived_n=int(arch_n))
    verdict["subsample_note"] = ("the 3,000-lead rows below apply 79_score_panels.subsample, which governs the GPU "
                                "panels (borzoi, chrombpnet); the archived alphagenome_atac row is all 11,896 leads, "
                                "so the reproduction guard compares on that set")
    # The written prediction is about the SAME variants, so it is evaluated on all served leads - the set the
    # archived row itself uses - with the 3,000-lead subsample reported beside it for reference only.
    # Evaluated on the variants BOTH models scored. The Atlas serves no indels, so its set is a subset of the
    # API's (10,966 of 11,896); comparing auROCs computed on the two different sets is the very error the
    # reproduction guard exists to catch.
    paired = both_scored(pd.to_numeric(df["atlas_atac"], errors="coerce").to_numpy(),
                         pd.to_numeric(df["ag_api_atac"], errors="coerce").to_numpy())
    yp = (pd.to_numeric(df["sign"], errors="coerce").to_numpy() > 0).astype(int)
    gp = df["chr"].astype(str).to_numpy()
    paired_atlas = score_one(yp, pd.to_numeric(df["atlas_atac"], errors="coerce").to_numpy(), gp, paired, "atlas_atac_paired")
    paired_api = score_one(yp, pd.to_numeric(df["ag_api_atac"], errors="coerce").to_numpy(), gp, paired, "ag_api_atac_paired")
    verdict["paired_same_variants"] = {"n": int(paired.sum()), "atlas_auroc": paired_atlas["auroc"],
                                       "api_auroc": paired_api["auroc"],
                                       "delta": float(abs(paired_atlas["auroc"] - paired_api["auroc"]))}
    delta_all = verdict["paired_same_variants"]["delta"]
    sub_atlas = verdict["atlas_atac_on_79_subsample_3000"]["auroc"]
    sub_api = verdict["ag_api_atac_on_79_subsample_3000"]["auroc"]
    verdict["prediction_P3a_caqtl"] = {
        "written": "|Atlas ATAC auROC - API ATAC auROC| < 0.02 on the same variants",
        "observed_delta_on_the_same_variants": float(delta_all),
        "n_paired": int(paired.sum()),
        "n_atlas_all": int(atlas_full["n"]), "n_api_all": int(api_full["n"]),
        "note_on_sets": "the Atlas serves no indels, so 930 API-scored leads have no Atlas score; the paired set is the comparison",
        "holds": bool(delta_all < 0.02),
        "observed_delta_on_79_subsample_reference_only": float(sub_atlas - sub_api) if sub_api == sub_api else None}
    json.dump(verdict, (TABLES / "external_benchmark_caqtl.json").open("w"), indent=1, default=float)
    la.log(f"caQTL: Atlas ATAC {atlas_full['auroc']:.4f} [{atlas_full['ci_lo']:.3f},{atlas_full['ci_hi']:.3f}] n={atlas_full['n']}; API recomputed {api_full['auroc']:.4f}; guard {verdict['G_api_reproduces_archive']['pass']} ({verdict['G_api_reproduces_archive']['reason']})")


def both_scored(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rows scored by BOTH models: the only set on which a paired auROC comparison is meaningful."""
    return (~np.isnan(np.asarray(a, float))) & (~np.isnan(np.asarray(b, float)))


def reproduction_guard(recomputed_auroc: float, recomputed_n: int, archived_auroc: float,
                       archived_n: int, tol: float = 0.005) -> dict:
    """Does this script's re-implementation reproduce the archived AlphaGenome-API caQTL number?

    A reproduction is only a reproduction on the same variant set, so the variant count must match the
    archive's own n. The archived alphagenome_atac row is all 11,896 served leads; the 3,000/4,000-lead
    subsample in 79_score_panels.subsample applies to the GPU panels (borzoi, chrombpnet), not to this row.
    """
    if recomputed_n != archived_n:
        return {"pass": False, "delta": None, "recomputed_n": int(recomputed_n), "archived_n": int(archived_n),
                "reason": f"different variant set: recomputed n {recomputed_n} != archived n {archived_n}"}
    delta = abs(float(recomputed_auroc) - float(archived_auroc))
    return {"pass": bool(delta < tol), "delta": delta, "recomputed_n": int(recomputed_n),
            "archived_n": int(archived_n), "tol": tol,
            "reason": "reproduced" if delta < tol else f"auROC differs by {delta:.6f} (tolerance {tol})"}


def gtex_intron_key(junction_start: int, junction_end: int) -> tuple:
    """Atlas half-open junction -> GTEx closed intron (shared with the joint_v2 join)."""
    return la.leafcutter_intron_key(junction_start, junction_end)


def atlas_junction_key(intron_start: int, intron_end: int) -> tuple:
    """GTEx (closed) intron coordinates -> Atlas (half-open) junction coordinates."""
    return (int(intron_start), int(intron_end) - 1)


def reconstruct_cluster(target: tuple, junction_keys: list) -> list:
    """Sibling junctions of `target`: those sharing its donor (start) or acceptor (end).

    GTEx's public sQTL release lists only the top phenotype per LeafCutter cluster (10,686 of 10,853 liver
    clusters are singletons in Liver.v8.sgenes, and all 1,414 clusters in the significant-pairs file are),
    so the cluster the slope is defined against is not distributed. LeafCutter's own definition - junctions
    sharing a splice site - is reconstructed here over the Atlas junction set for the same variant. It is a
    reconstruction, not GTEx's cluster, and every number derived from it carries that label.
    """
    ts, te = target
    return [(s, e) for (s, e) in junction_keys if (s, e) != target and (s == ts or e == te)]


def leafcutter_contrast(variant_uid: str, target: tuple, cluster: list, junc: dict) -> float:
    """Atlas junction quantile of the target intron minus the mean over the other introns in its cluster.

    GTEx sQTL slope is the ALT-allele effect on this intron's excision proportion WITHIN its LeafCutter
    cluster. An absolute junction quantile is a different estimand: an allele that lifts every junction in
    the cluster equally moves no proportion. Siblings the Atlas did not score are skipped rather than
    counted as zero, which would drag the contrast toward the target's own value.
    """
    t = junc.get((variant_uid, target[0], target[1]))
    if t is None or t != t:
        return math.nan
    others = [junc[(variant_uid, s, e)] for (s, e) in cluster
              if (s, e) != target and (variant_uid, s, e) in junc and junc[(variant_uid, s, e)] == junc[(variant_uid, s, e)]]
    if not others:
        return math.nan
    return float(t - np.mean(others))


def run_sqtl() -> None:
    import gzip
    sg = pd.read_csv(SQTL_DIR / "Liver.v8.sgenes.txt.gz", sep="\t", dtype=str)
    pairs = pd.read_csv(SQTL_DIR / "Liver.v8.sqtl_signifpairs.txt.gz", sep="\t", dtype=str)
    pairs["cluster"] = pairs["phenotype_id"].str.split(":").str[3]
    pairs["pval_nominal"] = pd.to_numeric(pairs["pval_nominal"], errors="coerce")
    # cluster membership comes from the FULL pair table; taking it after the head(10) subset would drop
    # sibling introns and make the within-cluster contrast a contrast against an arbitrary subset
    all_intr = pairs["phenotype_id"].str.split(":", expand=True)
    cluster_members = (pd.DataFrame({"cluster": pairs["cluster"],
                                     "s": all_intr[1].astype(int), "e": all_intr[2].astype(int)})
                       .drop_duplicates().groupby("cluster").apply(lambda g: list(zip(g["s"], g["e"])), include_groups=False).to_dict())
    pairs = pairs.sort_values("pval_nominal").groupby("cluster").head(10)
    def uid(v):
        c, p, r, a, _ = v.split("_")
        return f"{c}:{p}:{r}:{a}" if len(r) == 1 and len(a) == 1 else None
    pairs["variant_uid"] = pairs["variant_id"].map(uid)
    pairs = pairs.dropna(subset=["variant_uid"])
    pairs[["intron_start", "intron_end"]] = pairs["phenotype_id"].str.split(":", expand=True)[[1, 2]].astype(int)
    pairs["gene"] = pairs["phenotype_id"].str.split(":").str[4].str.split(".").str[0]
    # Atlas junction predictions per variant × junction × liver tracks
    junc = {}
    usage = {}
    for chunk in sorted((RAW / "atlas_sqtl").glob("chunk_*")):
        res = aa.load_archive(chunk)
        a = res.get("SPLICE_JUNCTIONS")
        if a is not None and a.shape[0] and "variant" in a.obs:
            obs = a.obs.reset_index(drop=True); mask = liver_track_mask(a.var.reset_index(drop=True), None)
            q = np.asarray(a.layers["quantiles"], np.float32)[:, mask]
            for i in range(len(obs)):
                u = aa.variant_uid_from_str(str(obs.at[i, "variant"]))
                js, je = int(float(obs.at[i, "junction_Start"])) if "junction_Start" in obs else -1, int(float(obs.at[i, "junction_End"])) if "junction_End" in obs else -1
                # stored in GTEx intron coordinates so cluster membership joins exactly
                junc[(u, *gtex_intron_key(js, je))] = float(np.nanmedian(q[i]))
        a = res.get("SPLICE_SITE_USAGE")
        if a is not None and a.shape[0] and "variant" in a.obs:
            obs = a.obs.reset_index(drop=True); mask = liver_track_mask(a.var.reset_index(drop=True), None)
            q = np.asarray(a.layers["quantiles"], np.float32)[:, mask]
            for i in range(len(obs)):
                u = aa.variant_uid_from_str(str(obs.at[i, "variant"])); g = str(obs.at[i, "gene_id"]).split(".")[0] if "gene_id" in obs else ""
                usage[(u, g)] = max(usage.get((u, g), 0.0), float(np.nanmax(np.abs(q[i]))))
    junc_by_uid = collections.defaultdict(list)
    for (u, s_, e_) in junc:
        junc_by_uid[u].append((s_, e_))
    rows = []
    for _, r in pairs.iterrows():
        jq = junc.get((r["variant_uid"], r["intron_start"], r["intron_end"]), math.nan)
        target = (r["intron_start"], r["intron_end"])
        sibs = reconstruct_cluster(target, junc_by_uid.get(r["variant_uid"], []))
        contrast = leafcutter_contrast(r["variant_uid"], target, [target] + sibs, junc)
        rows.append({"variant_uid": r["variant_uid"], "phenotype_id": r["phenotype_id"], "gene": r["gene"], "slope": float(r["slope"]), "pval_nominal": r["pval_nominal"],
                     "is_sgene_lead": r["variant_id"] in set(sg["variant_id"]), "atlas_junction_liver_quantile": jq,
                     "atlas_leafcutter_contrast": contrast, "n_reconstructed_siblings": len(sibs),
                     "atlas_usage_max_abs_quantile": usage.get((r["variant_uid"], r["gene"]), math.nan)})
    la.write_tsv_once(TABLES / "external_benchmark_sqtl.tsv", rows, list(rows[0].keys()))
    d = pd.DataFrame(rows); m = d.dropna(subset=["atlas_junction_liver_quantile"])
    m = m[m["atlas_junction_liver_quantile"] != 0]
    conc = float(((m["slope"] > 0) == (m["atlas_junction_liver_quantile"] > 0)).mean()) if len(m) else math.nan
    cm = d.dropna(subset=["atlas_leafcutter_contrast"])
    cm = cm[cm["atlas_leafcutter_contrast"] != 0]
    conc_contrast = float(((cm["slope"] > 0) == (cm["atlas_leafcutter_contrast"] > 0)).mean()) if len(cm) else math.nan
    cm_lead = cm[cm["is_sgene_lead"]]
    lead = m[m["is_sgene_lead"]]
    out = {"n_pairs": int(len(d)), "n_with_matching_junction": int(len(m)), "sign_concordance_matching_junction": conc,
           "sign_concordance_sgene_leads": float(((lead["slope"] > 0) == (lead["atlas_junction_liver_quantile"] > 0)).mean()) if len(lead) else math.nan, "n_sgene_leads_matched": int(len(lead)),
           "spearman_abs_slope_vs_abs_junction_quantile": float(spearmanr(m["slope"].abs(), m["atlas_junction_liver_quantile"].abs()).correlation) if len(m) > 10 else math.nan,
           "spearman_abs_slope_vs_usage": float(spearmanr(d.dropna(subset=["atlas_usage_max_abs_quantile"])["slope"].abs(), d.dropna(subset=["atlas_usage_max_abs_quantile"])["atlas_usage_max_abs_quantile"]).correlation) if d["atlas_usage_max_abs_quantile"].notna().sum() > 10 else math.nan,
           "n_with_leafcutter_contrast": int(len(cm)),
           "sign_concordance_leafcutter_contrast": conc_contrast,
           "sign_concordance_leafcutter_contrast_sgene_leads": float(((cm_lead["slope"] > 0) == (cm_lead["atlas_leafcutter_contrast"] > 0)).mean()) if len(cm_lead) else math.nan,
           "n_sgene_leads_contrast": int(len(cm_lead)),
           "spearman_slope_vs_contrast": float(spearmanr(cm["slope"], cm["atlas_leafcutter_contrast"]).correlation) if len(cm) > 10 else math.nan,
           "cluster_source": ("RECONSTRUCTED. GTEx's public sQTL release lists only the top phenotype per LeafCutter "
                              "cluster (all 1,414 clusters in the significant-pairs file are singletons; 10,686 of 10,853 "
                              "in sgenes), so the cluster the slope is defined against is not distributed. Siblings are "
                              "reconstructed as Atlas-scored junctions for the same variant sharing the target's donor or "
                              "acceptor, which is LeafCutter's own definition applied to a model-defined junction set."),
           "estimand": ("GTEx slope is a LeafCutter WITHIN-CLUSTER excision ratio; the absolute junction quantile is a "
                        "different quantity and its concordance is reported only as the mismatched comparison. The matched "
                        "quantity is atlas_leafcutter_contrast = target junction quantile - mean(sibling junction quantiles)."),
           "exposure": "documented_overlap: GTEx is an AlphaGenome training source; this is a ceiling/sanity number, not evidence",
           "junction_join": ("Atlas junction coordinates are half-open, GTEx LeafCutter introns are closed; junctions are "
                             "stored as (start, end+1) so the join is exact. The earlier +/-1 tolerance is removed."),
           "orientation": "GTEx slope is per ALT allele of the b38 variant id; Atlas quantile is ref->alt: same frame, no flip",
           "prediction_written": "sign concordance > 0.7"}
    json.dump(out, (TABLES / "external_benchmark_sqtl.json").open("w"), indent=1, default=float)
    la.log(f"sQTL: {out}")


def run_mpra() -> None:
    lib = pd.read_csv(TABLES / "mpra_library_intervals.tsv", sep="\t", dtype=str)
    chans = MPRA_CHANNELS
    rows = []
    for _, r in lib.iterrows():
        d = RAW / "atlas_mpra_intervals" / r["interval"].replace(":", "_")
        if not (d / "request.json").exists():
            continue
        res = aa.load_archive(d)
        centre = int(r["centre_pos_hg38"])
        row = {"interval": r["interval"], "centre_pos_hg38": centre, "allele_known": r["allele_known"], "ref_hg38": r["ref_hg38"], "alt_hg38": r["alt_hg38"],
               **{c: r[c] for c in lib.columns if c.startswith("dav_") or c.startswith("log2fc_")}}
        for label, (scorer, name_filter) in chans.items():
            a = res.get(scorer)
            if a is None or a.shape[0] == 0 or "variant" not in a.obs:
                row[f"{label}_centre_max_abs"] = math.nan; row[f"{label}_allele_signed"] = math.nan; continue
            obs = a.obs.reset_index(drop=True); var = a.var.reset_index(drop=True)
            if scorer == "AVI_SCORE":
                mask = np.ones(len(var), bool); vals = np.asarray(a.X, np.float32)
            else:
                if name_filter is None:
                    mask = var["ontology_curie"].isin(LIVER).values
                else:
                    mask = (var["biosample_name"].astype(str).str.contains(name_filter.split(".*")[0].split("|")[0], case=False).values
                            & (var["name"].astype(str).str.contains("H3K27ac").values if "H3K27ac" in name_filter else True))
                vals = np.asarray(a.layers["quantiles"], np.float32)
            if not mask.any():
                row[f"{label}_centre_max_abs"] = math.nan; row[f"{label}_allele_signed"] = math.nan; continue
            best, signed = math.nan, math.nan
            for i in range(len(obs)):
                u = aa.variant_uid_from_str(str(obs.at[i, "variant"]))
                if int(u.split(":")[1]) != centre:
                    continue
                v = float(np.nanmedian(vals[i][mask]))
                if best != best or abs(v) > best:
                    best = abs(v)
                if r["allele_known"] == "True" and u.split(":")[2] == r["ref_hg38"] and u.split(":")[3] == r["alt_hg38"]:
                    signed = v
            row[f"{label}_centre_max_abs"] = best; row[f"{label}_allele_signed"] = signed
        rows.append(row)
    la.write_tsv_once(TABLES / "external_benchmark_mpra_intervals.tsv", rows, list(rows[0].keys()))
    score_mpra(pd.DataFrame(rows), chans)


def rescore_mpra() -> None:
    """Score from the deposited per-interval table, without re-reading the 176 GB of interval archives."""
    d = pd.read_csv(TABLES / "external_benchmark_mpra_intervals.tsv", sep="\t")
    chans = MPRA_CHANNELS
    score_mpra(d, chans)


def magnitude_control_auc(stim_abs_effect: np.ndarray, const_abs_effect: np.ndarray) -> float:
    """Discrimination of stimulus-specific from constitutive DAVs using |MPRA effect| alone.

    Reported beside every stimulus-specific auROC. A value near 0.5 means the two groups have comparable
    measured effects, so a below-chance model auROC cannot be blamed on the model merely tracking effect
    size; a value below 0.5 means constitutive DAVs really do have larger effects and part of the model's
    below-chance result follows from that rather than from context-specific regulation.
    """
    a = np.asarray(stim_abs_effect, float); b = np.asarray(const_abs_effect, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if a.size == 0 or b.size == 0:
        return float("nan")
    y = np.concatenate([np.ones(a.size), np.zeros(b.size)])
    s = np.concatenate([a, b])
    return auroc(y, s, np.ones(y.size, bool))


def score_mpra(d: pd.DataFrame, chans: dict) -> None:
    signed = scorer_is_signed()
    d["chr"] = d["interval"].str.split(":").str[0]; d["block"] = d["chr"] + ":" + (d["centre_pos_hg38"] // 1_000_000).astype(str)
    out_rows = []
    for ctx in ("HepG2_ctrl", "HepG2_PAOA", "LX2_ctrl", "LX2_TGFB"):
        y = (d[f"dav_{ctx}"].astype(str) == "True").astype(int).to_numpy()
        for label in chans:
            s = pd.to_numeric(d[f"{label}_centre_max_abs"], errors="coerce").to_numpy(); ok = ~np.isnan(s)
            r = score_one(y, s, d["block"].to_numpy(), ok, label); r["context"] = ctx; r["n_pos"] = int(y[ok].sum()); out_rows.append(r)
        # stimulus-specific DAVs vs constitutive DAVs (same cell model): the static model should not separate these
        base = "HepG2_ctrl" if ctx.startswith("HepG2") else "LX2_ctrl"
        if ctx != base:
            stim_only = (d[f"dav_{ctx}"].astype(str) == "True") & (d[f"dav_{base}"].astype(str) != "True")
            both = (d[f"dav_{ctx}"].astype(str) == "True") & (d[f"dav_{base}"].astype(str) == "True")
            sub = d[stim_only | both]; y2 = stim_only[stim_only | both].astype(int).to_numpy()
            fc = pd.to_numeric(sub[f"log2fc_{ctx}"], errors="coerce").abs().to_numpy()
            mag = magnitude_control_auc(fc[y2 == 1], fc[y2 == 0])
            for label in ("hepg2_atac", "hepg2_h3k27ac", "avi"):
                s = pd.to_numeric(sub[f"{label}_centre_max_abs"], errors="coerce").to_numpy(); ok = ~np.isnan(s)
                r = score_one(y2, s, sub["block"].to_numpy(), ok, f"{label}|stimulus_specific_vs_constitutive")
                r["context"] = ctx; r["n_pos"] = int(y2[ok].sum())
                r["magnitude_control_auc"] = mag
                r["magnitude_explains"] = ("yes: constitutive DAVs have larger measured effects" if mag == mag and mag < 0.45 else
                                           "no: the two groups have comparable measured effects" if mag == mag and mag >= 0.45 else "unknown")
                out_rows.append(r)
        # allelic sign where alleles are known
        k = d[(d["allele_known"].astype(str) == "True") & (d[f"dav_{ctx}"].astype(str) == "True")]
        for label in chans:
            s = pd.to_numeric(k[f"{label}_allele_signed"], errors="coerce"); l2 = pd.to_numeric(k[f"log2fc_{ctx}"], errors="coerce"); ok = s.notna() & l2.notna() & (s != 0)
            conc = float(((s[ok] > 0) == (l2[ok] > 0)).mean()) if ok.sum() else math.nan
            out_rows.append({"context": ctx, "model": f"{label}|allelic_sign", "n": int(ok.sum()), "sign_concordance": conc,
                             "marginal_expected_concordance": la.marginal_expected_concordance(s[ok], l2[ok]) if ok.sum() else math.nan,
                             "sign_validity": sign_note(chans[label][0], signed)})
    la.write_tsv_once(TABLES / "external_benchmark_mpra.tsv", out_rows, sorted({k for r in out_rows for k in r}))
    json.dump({"n_intervals_scored": int(len(d)), "rule": "allele-agnostic centre-position max |quantile|; 1-Mb block bootstrap; label permutation",
               "prediction_written": "HepG2 auROC 0.60-0.70 with AVI ~ best single channel; stimulus-specific vs constitutive auROC within [0.45, 0.58]"},
              (TABLES / "external_benchmark_mpra.json").open("w"), indent=1)
    la.log(f"MPRA: {len(d)} intervals scored, {len(out_rows)} result rows")


if __name__ == "__main__":
    {"caqtl": run_caqtl, "sqtl": run_sqtl, "mpra": run_mpra, "mpra_rescore": rescore_mpra}[sys.argv[1]]()
