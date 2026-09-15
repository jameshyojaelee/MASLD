#!/usr/bin/env python3
"""C1: matched endogenous signed-allele-effect model comparison.

Executes scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md section 5 (C1 only).
Prespecification: PRESPECIFICATION.md in the output directory (written before any statistic).

Endpoint 1  Currin caQTL peak leads, label beta_alt (ALT-dosage FastQTL slope).
Endpoint 2  union allelic-imbalance sites, label mean_log2_alt_over_ref, weight 1/se.

Everything is scored on the matched set only: variants carrying a score from every compared model.
Resampling unit is the 1-Mb block; 10,000 draws, seed 20260914, paired across models in a draw.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import time

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
ATLAS = os.path.join(ROOT, "GWAS/finemapping/results/alphagenome_atlas")

CBP_ADULT = os.path.join(
    SF, "adult_liver_chrombpnet/hepatocyte_5fold_v2/currin_calibration/currin_scored_primary.tsv.gz")
CBP_281367 = os.path.join(
    SF, "adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/currin_calibration/"
        "currin_scored_primary.tsv.gz")
DIRFEAT = os.path.join(SF, "direction_features_caqtl.tsv")
P3A_CHUNKS = os.path.join(ATLAS, "p3a-benchmarks-20260909T190214Z/raw/atlas_caqtl")
ASE_281367 = os.path.join(ATLAS, "p3-ase-20260909T192926Z/tables/allelic_sites.tsv")
ASE_244832 = os.path.join(ATLAS, "p3-ase-20260909T192926Z/tables/gse244832_allelic_sites.tsv")
TRACK0 = os.path.join(ATLAS, "run-20260909T153939Z/tables/prediction_records")

BOOT_SEED = 20260914
N_BOOT = 10000

# adult-life-stage liver tracks only. The archived var table marks UBERON:0002107 DNase-seq and
# CL:0000182 DNase-seq embryonic; the three liver ATAC tracks are adult.
AG_ATAC_TRACKS = ["UBERON:0001114 ATAC-seq", "UBERON:0001115 ATAC-seq", "UBERON:0002107 ATAC-seq"]
AG_DNASE_TRACKS = ["UBERON:0001114 DNase-seq", "UBERON:0001115 DNase-seq"]
T0_ATAC_COLS = [f"raw|primary_liver|{t}" for t in AG_ATAC_TRACKS]
T0_DNASE_COLS = [f"raw|primary_liver|{t}" for t in AG_DNASE_TRACKS]

NT_RANK = {"A": 0, "C": 1, "G": 2, "T": 3}

MODELS_E1 = [
    "chrombpnet_adult_hep",
    "chrombpnet_adult_hep_gse281367",
    "chrombpnet_hepg2",
    "borzoi_atac",
    "borzoi_dnase",
    "alphagenome_atac_liver",
    "alphagenome_dnase_liver",
    "control_allele_identity",
    "control_position",
]
TIER_B_MODELS = [
    "chrombpnet_adult_hep",
    "chrombpnet_adult_hep_gse281367",
    "alphagenome_atac_liver",
    "alphagenome_dnase_liver",
    "control_allele_identity",
    "control_position",
]
MODELS_E2 = [
    "chrombpnet_adult_hep",
    "chrombpnet_adult_hep_gse281367",
    "alphagenome_atac_liver",
    "control_allele_identity",
    "control_position",
]

REJ = "reject_at_fwer_0.05"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256(path: str, max_bytes: int = 400_000_000) -> str:
    if os.path.getsize(path) > max_bytes:
        return "not_hashed_too_large"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def block_of(chrom, pos) -> str:
    c = str(chrom)
    c = c if c.startswith("chr") else f"chr{c}"
    return f"{c}~{int(pos) // 1_000_000}"


def key_of(chrom, pos, ref, alt) -> str:
    c = str(chrom)
    c = c[3:] if c.startswith("chr") else c
    return f"{c}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


# --------------------------------------------------------------------------- statistics
def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    da = np.sqrt((a * a).sum())
    db = np.sqrt((b * b).sum())
    if da <= 0 or db <= 0:
        return np.nan
    return float((a * b).sum() / (da * db))


def signed_spearman(score: np.ndarray, label: np.ndarray) -> float:
    if len(score) < 3:
        return np.nan
    return _pearson(stats.rankdata(score), stats.rankdata(label))


def weighted_spearman(score: np.ndarray, label: np.ndarray, w: np.ndarray) -> float:
    """Spearman with observation weights: rank both vectors, then weighted Pearson on ranks."""
    if len(score) < 3:
        return np.nan
    rs = stats.rankdata(score)
    rl = stats.rankdata(label)
    sw = w.sum()
    if sw <= 0:
        return np.nan
    ms = (w * rs).sum() / sw
    ml = (w * rl).sum() / sw
    cov = (w * (rs - ms) * (rl - ml)).sum() / sw
    vs = (w * (rs - ms) ** 2).sum() / sw
    vl = (w * (rl - ml) ** 2).sum() / sw
    if vs <= 0 or vl <= 0:
        return np.nan
    return float(cov / np.sqrt(vs * vl))


def sign_concordance(score: np.ndarray, label: np.ndarray, w: np.ndarray | None = None) -> float:
    ok = (score != 0) & (label != 0) & np.isfinite(score) & np.isfinite(label)
    if ok.sum() == 0:
        return np.nan
    agree = (np.sign(score[ok]) == np.sign(label[ok])).astype(float)
    if w is None:
        return float(agree.mean())
    ww = w[ok]
    return float((agree * ww).sum() / ww.sum())


class BlockBootstrap:
    """Shared block-resampling draws so every model is evaluated on the same resampled rows."""

    def __init__(self, blocks: np.ndarray, n_boot: int, seed: int):
        order = np.argsort(blocks, kind="mergesort")
        if not np.array_equal(order, np.arange(len(blocks))):
            raise ValueError("rows must be pre-sorted by block")
        uniq, starts, sizes = np.unique(blocks, return_index=True, return_counts=True)
        self.starts = starts.astype(np.int64)
        self.sizes = sizes.astype(np.int64)
        self.n_blocks = len(uniq)
        self.n_boot = n_boot
        self.rng = np.random.default_rng(seed)

    def draws(self):
        for _ in range(self.n_boot):
            pick = self.rng.integers(0, self.n_blocks, size=self.n_blocks)
            s = self.sizes[pick]
            st = self.starts[pick]
            total = int(s.sum())
            ends = np.cumsum(s)
            offs = np.arange(total) - np.repeat(ends - s, s)
            yield np.repeat(st, s) + offs


def ci(vals: np.ndarray) -> tuple[float, float]:
    v = vals[np.isfinite(vals)]
    if len(v) == 0:
        return (np.nan, np.nan)
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))


def two_sided_boot_p(diffs: np.ndarray) -> tuple[float, str]:
    d = diffs[np.isfinite(diffs)]
    if len(d) == 0:
        return (np.nan, "no_valid_draws")
    p = 2.0 * min(float((d <= 0).mean()), float((d >= 0).mean()))
    p = min(1.0, p)
    floor = 1.0 / len(d)
    if p < floor:
        return (floor, f"floored_at_1/{len(d)}_draws")
    return (float(p), "")


# --------------------------------------------------------------------------- loaders
def load_cbp(path: str, tag: str) -> pd.DataFrame:
    cols = ["variant_id_hg38", "chr_x", "pos_hg38", "ref", "alt", "label", "beta_alt",
            "q_value", "p_nominal", "peak_id", "reference_match", "heldout_fold", "logfc",
            "jsd", "active_allele_quantile", "model_peak_overlap"]
    d = pd.read_csv(path, sep="\t", usecols=cols)
    d["key"] = [key_of(c, p, r, a) for c, p, r, a in zip(d.chr_x, d.pos_hg38, d.ref, d.alt)]
    d = d.rename(columns={"logfc": tag})
    log(f"  {tag}: {len(d)} scored SNVs, {(d.label == 1).sum()} with beta_alt, "
        f"reference_match false = {int((~d.reference_match.astype(bool)).sum())}")
    return d


def load_atlas_chunks() -> pd.DataFrame:
    import anndata as ad
    frames = []
    chunks = sorted(glob.glob(os.path.join(P3A_CHUNKS, "chunk_*")))
    log(f"  reading {len(chunks)} archived Atlas caQTL chunks")
    for d in chunks:
        per = {}
        for stem, tracks in (("ATAC", AG_ATAC_TRACKS), ("DNASE", AG_DNASE_TRACKS)):
            a = ad.read_h5ad(os.path.join(d, f"{stem}.h5ad"))
            names = a.var["name"].astype(str).values
            idx = [int(np.where(names == t)[0][0]) for t in tracks]
            X = np.asarray(a.X, dtype=float)[:, idx]
            Q = np.asarray(a.layers["quantiles"], dtype=float)[:, idx]
            v = pd.Series(a.obs["variant"].astype(str).values)
            parts = v.str.replace(">", ":", regex=False).str.split(":", expand=True)
            per[stem] = pd.DataFrame({
                "key": [key_of(c, p, r, a2) for c, p, r, a2 in
                        zip(parts[0], parts[1], parts[2], parts[3])],
                f"{stem}_raw": np.nanmean(X, axis=1),
                f"{stem}_q": np.nanmean(Q, axis=1)})
        m = per["ATAC"].merge(per["DNASE"], on="key", how="outer", validate="one_to_one")
        frames.append(m)
    out = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["key"])
    out = out.rename(columns={"ATAC_raw": "alphagenome_atac_liver",
                              "DNASE_raw": "alphagenome_dnase_liver",
                              "ATAC_q": "alphagenome_atac_liver_quantile",
                              "DNASE_q": "alphagenome_dnase_liver_quantile"})
    log(f"  archived Atlas caQTL scores: {len(out)} unique variants")
    return out


def load_track0_ase(uids: list[str]) -> pd.DataFrame:
    import pyarrow.parquet as pq
    want = set(uids)
    frames = {}
    for stem, cols in (("ATAC", T0_ATAC_COLS), ("DNASE", T0_DNASE_COLS)):
        f = pq.ParquetFile(os.path.join(TRACK0, f"{stem}.parquet"))
        keep = []
        for batch in f.iter_batches(batch_size=50_000, columns=["variant_uid"] + cols):
            b = batch.to_pandas()
            keep.append(b[b.variant_uid.isin(want)])
        d = pd.concat(keep, ignore_index=True).drop_duplicates(subset=["variant_uid"])
        d[f"alphagenome_{stem.lower()}_liver"] = d[cols].mean(axis=1)
        frames[stem] = d[["variant_uid", f"alphagenome_{stem.lower()}_liver"]]
        log(f"  Track0 {stem}: {len(d)} of {len(want)} ASE sites served")
    out = frames["ATAC"].merge(frames["DNASE"], on="variant_uid", how="outer")
    return out


# --------------------------------------------------------------------------- endpoint runner
def run_endpoint(df: pd.DataFrame, models: list[str], label_col: str, weight_col: str | None,
                 name: str, out_rows: list, boot_store: dict) -> pd.DataFrame:
    df = df.sort_values("block_1mb", kind="mergesort").reset_index(drop=True)
    blocks = df["block_1mb"].values
    lab = df[label_col].values.astype(float)
    w = df[weight_col].values.astype(float) if weight_col else None
    n_blocks = int(pd.unique(blocks).size)
    log(f"[{name}] n={len(df)} variants, {n_blocks} 1-Mb blocks, {len(models)} models")

    point = {}
    for m in models:
        s = df[m].values.astype(float)
        rho = weighted_spearman(s, lab, w) if w is not None else signed_spearman(s, lab)
        point[m] = rho
        out_rows.append(dict(endpoint=name, model=m, metric="signed_spearman", value=rho,
                             n_variants=len(df), n_blocks=n_blocks))
        out_rows.append(dict(endpoint=name, model=m, metric="sign_concordance",
                             value=sign_concordance(s, lab, w),
                             n_variants=len(df), n_blocks=n_blocks))
        if w is not None:
            out_rows.append(dict(endpoint=name, model=m, metric="signed_spearman_unweighted",
                                 value=signed_spearman(s, lab),
                                 n_variants=len(df), n_blocks=n_blocks))
            out_rows.append(dict(endpoint=name, model=m, metric="sign_concordance_unweighted",
                                 value=sign_concordance(s, lab),
                                 n_variants=len(df), n_blocks=n_blocks))

    log(f"[{name}] point: " + ", ".join(
        f"{m}={point[m]:.4f}" for m in models if np.isfinite(point[m])))

    draws = {m: np.full(N_BOOT, np.nan) for m in models}
    scores = {m: df[m].values.astype(float) for m in models}
    bb = BlockBootstrap(blocks, N_BOOT, BOOT_SEED)
    t0 = time.time()
    for i, idx in enumerate(bb.draws()):
        y = lab[idx]
        ww = w[idx] if w is not None else None
        for m in models:
            s = scores[m][idx]
            draws[m][i] = weighted_spearman(s, y, ww) if ww is not None else signed_spearman(s, y)
        if (i + 1) % 2000 == 0:
            log(f"  [{name}] bootstrap {i+1}/{N_BOOT} ({time.time()-t0:.0f}s)")
    boot_store[name] = draws

    for m in models:
        lo, hi = ci(draws[m])
        out_rows.append(dict(endpoint=name, model=m, metric="signed_spearman_boot_ci_lo",
                             value=lo, n_variants=len(df), n_blocks=n_blocks))
        out_rows.append(dict(endpoint=name, model=m, metric="signed_spearman_boot_ci_hi",
                             value=hi, n_variants=len(df), n_blocks=n_blocks))
    return df


def contrasts(draws: dict, pairs: list[tuple[str, str]], name: str, rows: list,
              point_lookup: dict) -> None:
    for a, b in pairs:
        if a not in draws or b not in draws:
            continue
        d = draws[a] - draws[b]
        lo, hi = ci(d)
        p, note = two_sided_boot_p(d)
        rows.append(dict(endpoint=name, contrast=f"{a} - {b}",
                         point=point_lookup[(name, a)] - point_lookup[(name, b)],
                         boot_mean=float(np.nanmean(d)), ci_lo=lo, ci_hi=hi,
                         excludes_zero=bool(lo > 0 or hi < 0),
                         boot_p_two_sided=p, note=note,
                         n_draws=int(np.isfinite(d).sum())))


# --------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out = args.outdir
    tables = os.path.join(out, "tables")
    os.makedirs(tables, exist_ok=True)

    log("=== loading model scores ===")
    cbp_a = load_cbp(CBP_ADULT, "chrombpnet_adult_hep")
    cbp_b = load_cbp(CBP_281367, "chrombpnet_adult_hep_gse281367")

    dirfeat = pd.read_csv(DIRFEAT, sep="\t")
    dirfeat["key"] = dirfeat["canon"].astype(str)
    dirfeat = dirfeat.rename(columns={
        "cbp_logfc_labelframe": "chrombpnet_hepg2",
        "borzoi_atac_delta_labelframe": "borzoi_atac",
        "borzoi_dnase_delta_labelframe": "borzoi_dnase",
        "ag_atac_delta_labelframe": "ag_atac_src79",
        "ag_dnase_delta_labelframe": "ag_dnase_src79"})
    log(f"  direction_features_caqtl: {len(dirfeat)} peak leads; "
        f"hepg2={int(dirfeat.chrombpnet_hepg2.notna().sum())}, "
        f"borzoi_atac={int(dirfeat.borzoi_atac.notna().sum())}, "
        f"ag_atac={int(dirfeat.ag_atac_src79.notna().sum())}")

    atlas = load_atlas_chunks()

    # ---- endpoint 1 assembly -------------------------------------------------
    lab1 = cbp_a[cbp_a.label == 1].copy()
    lab1 = lab1[["key", "chr_x", "pos_hg38", "ref", "alt", "beta_alt", "q_value", "peak_id",
                 "heldout_fold", "chrombpnet_adult_hep", "model_peak_overlap"]]
    lab1 = lab1.merge(cbp_b[["key", "chrombpnet_adult_hep_gse281367"]], on="key", how="left")
    lab1 = lab1.merge(dirfeat[["key", "chrombpnet_hepg2", "borzoi_atac", "borzoi_dnase",
                               "ag_atac_src79", "ag_dnase_src79", "effect", "sign",
                               "feature_id", "maf", "tss_or_peak_dist"]], on="key", how="left")
    lab1 = lab1.merge(atlas, on="key", how="left")
    lab1["block_1mb"] = [block_of(c, p) for c, p in zip(lab1.chr_x, lab1.pos_hg38)]
    lab1["control_allele_identity"] = [
        NT_RANK.get(str(a).upper(), np.nan) - NT_RANK.get(str(r).upper(), np.nan)
        for r, a in zip(lab1.ref, lab1.alt)]
    lab1["control_position"] = (lab1.pos_hg38 % 1000) / 1000.0 - 0.5

    # ---- orientation and identity QC ----------------------------------------
    qc = {"n_currin_leads_with_beta_alt": int(len(lab1))}
    shared = lab1[lab1.effect.notna()]
    qc["label_vs_direction_features_n"] = int(len(shared))
    qc["label_vs_direction_features_sign_agreement"] = float(
        (np.sign(shared.beta_alt) == np.sign(shared.effect)).mean())
    qc["label_vs_direction_features_spearman"] = signed_spearman(
        shared.beta_alt.values.astype(float), shared.effect.values.astype(float))
    qc["label_vs_direction_features_same_peak_share"] = float(
        (shared.peak_id == shared.feature_id).mean())
    ags = lab1[lab1.ag_atac_src79.notna() & lab1.alphagenome_atac_liver.notna()]
    qc["atlas_archive_vs_src79_ag_atac_n"] = int(len(ags))
    qc["atlas_archive_vs_src79_ag_atac_spearman"] = signed_spearman(
        ags.alphagenome_atac_liver.values.astype(float), ags.ag_atac_src79.values.astype(float))
    qc["atlas_archive_vs_src79_ag_atac_sign_agreement"] = float(
        (np.sign(ags.alphagenome_atac_liver) == np.sign(ags.ag_atac_src79)).mean())
    cb = lab1[lab1.chrombpnet_adult_hep.notna() & lab1.chrombpnet_adult_hep_gse281367.notna()]
    qc["chrombpnet_family_agreement_spearman"] = signed_spearman(
        cb.chrombpnet_adult_hep.values.astype(float),
        cb.chrombpnet_adult_hep_gse281367.values.astype(float))
    log("  QC " + json.dumps(qc, indent=1))

    tierA = lab1.dropna(subset=MODELS_E1 + ["beta_alt"]).copy()
    tierB = lab1.dropna(subset=TIER_B_MODELS + ["beta_alt"]).copy()
    log(f"  tier A matched n={len(tierA)} blocks={tierA.block_1mb.nunique()}")
    log(f"  tier B matched n={len(tierB)} blocks={tierB.block_1mb.nunique()}")

    # ---- endpoint 2 assembly -------------------------------------------------
    a1 = pd.read_csv(ASE_281367, sep="\t")
    a2 = pd.read_csv(ASE_244832, sep="\t")
    a1["dataset"] = "GSE281367"
    a2["dataset"] = "GSE244832"
    both = pd.concat([a1, a2], ignore_index=True)
    rec = []
    for uid, g in both.groupby("uid"):
        prec = 1.0 / g.se_log2 ** 2
        wsum = float(prec.sum())
        rec.append(dict(uid=uid, chrom=g.chrom.iloc[0], pos=int(g.pos.iloc[0]),
                        ref=g.ref.iloc[0], alt=g.alt.iloc[0],
                        mean_log2_alt_over_ref=float((g.mean_log2_alt_over_ref * prec).sum() / wsum),
                        se_log2=float(np.sqrt(1.0 / wsum)),
                        n_datasets=int(len(g)), datasets=",".join(sorted(g.dataset)),
                        n_het_donors=int(g.n_het_donors.sum()),
                        total_reads=int(g.total_reads.sum())))
    ase = pd.DataFrame(rec)
    ase["key"] = [key_of(c, p, r, a) for c, p, r, a in zip(ase.chrom, ase.pos, ase.ref, ase.alt)]
    ase["block_1mb"] = [block_of(c, p) for c, p in zip(ase.chrom, ase.pos)]
    ase["precision_weight"] = 1.0 / ase.se_log2
    log(f"  ASE union: {len(ase)} sites, {ase.block_1mb.nunique()} blocks, "
        f"{int((ase.n_datasets == 2).sum())} in both datasets")

    t0tab = load_track0_ase(ase.uid.tolist())
    ase = ase.merge(t0tab.rename(columns={"variant_uid": "uid"}), on="uid", how="left")
    ase = ase.merge(cbp_a[["key", "chrombpnet_adult_hep", "label", "beta_alt"]]
                    .rename(columns={"label": "in_currin_positive_set",
                                     "beta_alt": "currin_beta_alt"}), on="key", how="left")
    ase = ase.merge(cbp_b[["key", "chrombpnet_adult_hep_gse281367"]], on="key", how="left")
    ase = ase.merge(dirfeat[["key", "chrombpnet_hepg2", "borzoi_atac", "borzoi_dnase"]],
                    on="key", how="left")
    ase["control_allele_identity"] = [
        NT_RANK.get(str(a).upper(), np.nan) - NT_RANK.get(str(r).upper(), np.nan)
        for r, a in zip(ase.ref, ase.alt)]
    ase["control_position"] = (ase.pos % 1000) / 1000.0 - 0.5
    cov_models = ["alphagenome_atac_liver", "alphagenome_dnase_liver", "chrombpnet_adult_hep",
                  "chrombpnet_adult_hep_gse281367", "chrombpnet_hepg2", "borzoi_atac",
                  "borzoi_dnase"]
    cov = {m: int(ase[m].notna().sum()) for m in cov_models}
    log(f"  ASE coverage per model: {cov}")
    ase.to_csv(os.path.join(tables, "endpoint2_ase_union_scores.tsv"), sep="\t", index=False)

    ase_m = ase.dropna(subset=MODELS_E2 + ["mean_log2_alt_over_ref", "se_log2"]).copy()
    log(f"  endpoint 2 matched n={len(ase_m)} blocks={ase_m.block_1mb.nunique()}")

    # ---- run -----------------------------------------------------------------
    rows, boot = [], {}
    tierA = run_endpoint(tierA, MODELS_E1, "beta_alt", None, "endpoint1_currin_tierA", rows, boot)
    tierB = run_endpoint(tierB, TIER_B_MODELS, "beta_alt", None, "endpoint1_currin_tierB",
                         rows, boot)
    ase_m = run_endpoint(ase_m, MODELS_E2, "mean_log2_alt_over_ref", "precision_weight",
                         "endpoint2_allelic_imbalance", rows, boot)

    tierA.to_csv(os.path.join(tables, "endpoint1_matched_tierA.tsv.gz"), sep="\t",
                 index=False, compression="gzip")
    tierB.to_csv(os.path.join(tables, "endpoint1_matched_tierB.tsv.gz"), sep="\t",
                 index=False, compression="gzip")
    ase_m.to_csv(os.path.join(tables, "endpoint2_matched.tsv"), sep="\t", index=False)

    metrics = pd.DataFrame(rows)
    metrics.to_csv(os.path.join(tables, "c1_model_metrics.tsv"), sep="\t", index=False)

    plook = {(r.endpoint, r.model): r.value for r in
             metrics[metrics.metric == "signed_spearman"].itertuples()}
    crows = []
    e1_pairs = [("chrombpnet_adult_hep", "alphagenome_atac_liver"),
                ("chrombpnet_adult_hep", "alphagenome_dnase_liver"),
                ("chrombpnet_adult_hep", "chrombpnet_hepg2"),
                ("chrombpnet_adult_hep", "borzoi_atac"),
                ("chrombpnet_adult_hep", "borzoi_dnase"),
                ("chrombpnet_adult_hep", "control_allele_identity"),
                ("alphagenome_atac_liver", "control_allele_identity"),
                ("borzoi_atac", "alphagenome_atac_liver"),
                ("chrombpnet_adult_hep_gse281367", "alphagenome_atac_liver")]
    contrasts(boot["endpoint1_currin_tierA"], e1_pairs, "endpoint1_currin_tierA", crows, plook)
    contrasts(boot["endpoint1_currin_tierB"],
              [p for p in e1_pairs if p[0] in TIER_B_MODELS and p[1] in TIER_B_MODELS],
              "endpoint1_currin_tierB", crows, plook)
    e2_pairs = [("chrombpnet_adult_hep", "alphagenome_atac_liver"),
                ("chrombpnet_adult_hep", "control_allele_identity"),
                ("alphagenome_atac_liver", "control_allele_identity")]
    contrasts(boot["endpoint2_allelic_imbalance"], e2_pairs, "endpoint2_allelic_imbalance",
              crows, plook)
    con = pd.DataFrame(crows)
    con.to_csv(os.path.join(tables, "c1_contrasts.tsv"), sep="\t", index=False)
    log("contrasts\n" + con.to_string())

    # ---- Holm over the confirmatory family of three --------------------------
    p1 = float(con[(con.endpoint == "endpoint1_currin_tierA")
                   & (con.contrast == "chrombpnet_adult_hep - alphagenome_atac_liver")
                   ].boot_p_two_sided.iloc[0])
    p2 = float(con[(con.endpoint == "endpoint2_allelic_imbalance")
                   & (con.contrast == "chrombpnet_adult_hep - alphagenome_atac_liver")
                   ].boot_p_two_sided.iloc[0])
    fam = [("C1_endpoint1_chrombpnet_vs_atlas_atac", p1),
           ("C1_endpoint2_chrombpnet_vs_atlas_atac", p2)]
    fam.sort(key=lambda t: t[1])
    m_family = 3
    holm, prev = [], 0.0
    for i, (k, p) in enumerate(fam):
        adj = max(prev, min(1.0, p * (m_family - i)))
        prev = adj
        holm.append({"member": k, "raw_p": p, "holm_step_threshold": 0.05 / (m_family - i),
                     "holm_adjusted_p": adj, REJ: bool(adj <= 0.05)})
    holm.append({"member": "E2_transfer_contrast", "raw_p": np.nan,
                 "holm_step_threshold": np.nan, "holm_adjusted_p": np.nan, REJ: False})
    hd = pd.DataFrame(holm)
    hd.to_csv(os.path.join(tables, "c1_holm_family.tsv"), sep="\t", index=False)
    log("holm\n" + hd.to_string())

    # ---- counts --------------------------------------------------------------
    counts = []
    for nm, d, ms in (("endpoint1_currin_tierA", tierA, MODELS_E1),
                      ("endpoint1_currin_tierB", tierB, TIER_B_MODELS),
                      ("endpoint2_allelic_imbalance", ase_m, MODELS_E2)):
        counts.append(dict(matched_set=nm, n_variants=len(d),
                           n_blocks_1mb=int(d.block_1mb.nunique()), models=",".join(ms)))
    counts.append(dict(matched_set="endpoint1_label_universe_currin_leads_with_beta",
                       n_variants=int(len(lab1)),
                       n_blocks_1mb=int(lab1.block_1mb.nunique()), models="label_only"))
    counts.append(dict(matched_set="endpoint2_ase_union_sites", n_variants=int(len(ase)),
                       n_blocks_1mb=int(ase.block_1mb.nunique()), models="label_only"))
    for m in ["chrombpnet_adult_hep", "chrombpnet_adult_hep_gse281367", "chrombpnet_hepg2",
              "borzoi_atac", "borzoi_dnase", "alphagenome_atac_liver", "alphagenome_dnase_liver"]:
        counts.append(dict(matched_set=f"endpoint1_coverage:{m}",
                           n_variants=int(lab1[m].notna().sum()),
                           n_blocks_1mb=int(lab1.loc[lab1[m].notna(), "block_1mb"].nunique()),
                           models=m))
    for m, n in cov.items():
        counts.append(dict(matched_set=f"endpoint2_coverage:{m}", n_variants=n,
                           n_blocks_1mb=int(ase.loc[ase[m].notna(), "block_1mb"].nunique()),
                           models=m))
    pd.DataFrame(counts).to_csv(os.path.join(tables, "c1_matched_counts.tsv"),
                                sep="\t", index=False)

    # ---- sensitivity: AlphaGenome quantile layer -----------------------------
    sens = []
    for nm, d in (("endpoint1_currin_tierA", tierA), ("endpoint1_currin_tierB", tierB)):
        for col in ("alphagenome_atac_liver_quantile", "alphagenome_dnase_liver_quantile"):
            sens.append(dict(endpoint=nm, model=col, metric="signed_spearman",
                             value=signed_spearman(d[col].values.astype(float),
                                                   d.beta_alt.values.astype(float)),
                             n_variants=len(d)))
        for col in ("ag_atac_src79", "ag_dnase_src79"):
            v = d[d[col].notna()]
            sens.append(dict(endpoint=nm, model=col, metric="signed_spearman",
                             value=signed_spearman(v[col].values.astype(float),
                                                   v.beta_alt.values.astype(float)),
                             n_variants=len(v)))
    pd.DataFrame(sens).to_csv(os.path.join(tables, "c1_alphagenome_sensitivity.tsv"),
                              sep="\t", index=False)

    with open(os.path.join(out, "qc_orientation.json"), "w") as fh:
        json.dump(qc, fh, indent=2)

    # ---- manifest ------------------------------------------------------------
    man = []
    for p in (CBP_ADULT, CBP_281367, DIRFEAT, ASE_281367, ASE_244832,
              os.path.join(TRACK0, "ATAC.parquet"), os.path.join(TRACK0, "DNASE.parquet"),
              os.path.join(P3A_CHUNKS, "chunk_00000", "ATAC.h5ad"),
              os.path.join(SF, "adult_liver_chrombpnet/hepatocyte_5fold_v2/currin_calibration/"
                               "currin_external_calibration.json")):
        man.append(dict(path=p, bytes=os.path.getsize(p), sha256=sha256(p)))
    man.append(dict(path=P3A_CHUNKS, bytes=-1,
                    sha256="directory_%d_chunks" % len(
                        glob.glob(os.path.join(P3A_CHUNKS, "chunk_*")))))
    pd.DataFrame(man).to_csv(os.path.join(out, "MANIFEST.tsv"), sep="\t", index=False)
    log("done")


if __name__ == "__main__":
    main()
