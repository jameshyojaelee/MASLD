#!/usr/bin/env python
"""A1. Reporter measurement reliability (GSE281364 MPRA).

Executes section 3 of scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md.

Quantities per context, per element:
  activity   a = log2((RNA+0.5)/(DNA+0.5)), element-level = mean over the ref and alt construct
  allele     d = a_alt - a_ref within a sample
  treatment  j = d_treated - d_control per cell line

Reliability = Pearson / Spearman between the two half-means of a balanced 2-versus-2
replicate partition, averaged over the three partitions {12|34, 13|24, 14|23}.
Spearman-Brown projection to four replicates is reported beside, never instead of, the
raw two-versus-two value.

Uncertainty: 10,000 bootstrap draws over the 239 Borzoi long-range blocks, seed 20260914.
The statistic is regenerated end to end inside every draw through the same function that
produced the observed value.

Usage: python a1_mpra_reliability.py <output_dir>
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import subprocess
import sys
import time
from multiprocessing import Pool

import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

PROJ = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BENCH = f"{PROJ}/Analysis/MASLD_Model_Benchmark"
VAL = f"{BENCH}/executions/gse281364-validated-reconstruction-21066470/validated"
SPLIT = f"{BENCH}/executions/gse281364-outcome-blind-splits-21068885/split"
BLOCKMAP = f"{BENCH}/executions/gse281364-borzoi-native-fixture-21083008/fixture/source_group_map.tsv"
FIXTURE = f"{BENCH}/executions/gse281364-dna-lm-common-fixture-21069076/fixture/sequence_manifest.tsv"

PC = 0.5
REPS = ["1", "2", "3", "4"]
SPLITS = [(("1", "2"), ("3", "4")), (("1", "3"), ("2", "4")), (("1", "4"), ("2", "3"))]
CONTEXTS = ["HepG2_control", "HepG2_PAOA", "LX2_control", "LX2_TGFb"]
CELL_PAIRS = {"HepG2": ("HepG2_PAOA", "HepG2_control"), "LX2": ("LX2_TGFb", "LX2_control")}
BOOT_SEED = 20260914
N_BOOT = int(os.environ.get("A1_N_BOOT", "10000"))
PERM_SEED = 20260914
N_PERM = int(os.environ.get("A1_N_PERM", "10000"))

INPUTS = [
    f"{VAL}/replicate_outcomes.tsv.gz",
    f"{VAL}/elements.tsv",
    f"{SPLIT}/groups.tsv",
    f"{SPLIT}/elements.tsv",
    BLOCKMAP,
    FIXTURE,
]


# ---------------------------------------------------------------- loading


def load():
    el = pd.read_csv(f"{VAL}/elements.tsv", sep="\t", dtype=str,
                     keep_default_na=False, na_values=[])
    ro = pd.read_csv(f"{VAL}/replicate_outcomes.tsv.gz", sep="\t", dtype=str,
                     keep_default_na=False, na_values=[])
    paired = el[el.pair_state == "paired_snv"].copy()
    keep = set(paired.element_id)
    grp = paired.set_index("element_id")["outer_locus_sequence_group_id"]

    bm = pd.read_csv(BLOCKMAP, sep="\t", dtype=str)
    g2b = bm.set_index("source_locus_group_id")["borzoi_long_range_group_id"]

    sub = ro[ro.element_id.isin(keep)].copy()
    qc = sub[sub.assay_state != "observed"].copy()
    o = sub[sub.assay_state == "observed"].copy()
    for c in ["DNA", "RNA", "n_barcodes"]:
        o[c] = pd.to_numeric(o[c])
    o["act"] = np.log2((o.RNA + PC) / (o.DNA + PC))
    return el, ro, paired, grp, g2b, o, qc


def build_matrices(o, elements):
    """Return dict keyed by (context, allele) -> DataFrame [elements x 4 replicates]."""
    out = {}
    for (ctx, allele), g in o.groupby(["context_id", "allele"]):
        for val in ["act", "DNA", "RNA", "n_barcodes"]:
            p = g.pivot_table(index="element_id", columns="experimental_replicate",
                              values=val, aggfunc="first")
            for r in REPS:
                if r not in p.columns:
                    p[r] = np.nan
            out[(ctx, allele, val)] = p[REPS].reindex(elements)
    return out


# ---------------------------------------------------------------- statistic


def half_means(M, partition_index):
    """M: (n,4) array of per-replicate values. Returns (n,2) half-means."""
    a, b = SPLITS[partition_index]
    ia = [REPS.index(x) for x in a]
    ib = [REPS.index(x) for x in b]
    return np.column_stack([M[:, ia].mean(axis=1), M[:, ib].mean(axis=1)])


def _pearson(x, y):
    xm = x - x.mean()
    ym = y - y.mean()
    dx = np.sqrt((xm * xm).sum())
    dy = np.sqrt((ym * ym).sum())
    if dx == 0 or dy == 0:
        return np.nan
    return float((xm * ym).sum() / (dx * dy))


def _spearman(x, y):
    return _pearson(stats.rankdata(x), stats.rankdata(y))


def sb_project(r, k_from=2, k_to=4):
    """Spearman-Brown. r is the correlation between two means of k_from replicates."""
    if not np.isfinite(r) or r >= 1.0:
        return np.nan
    r1 = r / (k_from - (k_from - 1) * r)
    if r1 <= -1 / (k_to - 1):
        return np.nan
    return float(k_to * r1 / (1 + (k_to - 1) * r1))


def reliability_from_halves(H, idx=None):
    """H: (n, 3, 2) half-means for the three partitions. Returns dict of statistics."""
    if idx is not None:
        H = H[idx]
    pear = np.array([_pearson(H[:, p, 0], H[:, p, 1]) for p in range(H.shape[1])])
    spear = np.array([_spearman(H[:, p, 0], H[:, p, 1]) for p in range(H.shape[1])])
    mp = float(np.nanmean(pear))
    ms = float(np.nanmean(spear))
    return {"pearson": mp, "spearman": ms,
            "projected_4rep": sb_project(mp),
            "projected_4rep_spearman": sb_project(ms),
            "pearson_by_partition": pear, "spearman_by_partition": spear}


# ---------------------------------------------------------------- bootstrap


def block_index(block_codes, n_blocks):
    order = np.argsort(block_codes, kind="stable")
    sorted_codes = block_codes[order]
    starts = np.searchsorted(sorted_codes, np.arange(n_blocks), side="left")
    ends = np.searchsorted(sorted_codes, np.arange(n_blocks), side="right")
    return order, starts, ends


def bootstrap_ci(H, block_codes, n_blocks, n_boot=N_BOOT, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    order, starts, ends = block_index(block_codes, n_blocks)
    draws = rng.integers(0, n_blocks, size=(n_boot, n_blocks))
    res = {k: np.full(n_boot, np.nan) for k in
           ["pearson", "spearman", "projected_4rep", "projected_4rep_spearman"]}
    for b in range(n_boot):
        sel = draws[b]
        idx = np.concatenate([order[starts[s]:ends[s]] for s in sel])
        st = reliability_from_halves(H, idx)
        for k in res:
            res[k][b] = st[k]
    out = {}
    for k, v in res.items():
        v = v[np.isfinite(v)]
        if len(v) == 0:
            out[f"{k}_lo"] = np.nan
            out[f"{k}_hi"] = np.nan
        else:
            out[f"{k}_lo"] = float(np.percentile(v, 2.5))
            out[f"{k}_hi"] = float(np.percentile(v, 97.5))
        out[f"{k}_n_finite_draws"] = int(len(v))
    return out


def _boot_worker(args):
    key, H, codes, nb = args
    t0 = time.time()
    ci = bootstrap_ci(H, codes, nb)
    ci["_boot_seconds"] = round(time.time() - t0, 1)
    return key, ci


# ---------------------------------------------------------------- main


def sha256(path, limit_gb=2):
    if os.path.getsize(path) > limit_gb * 1024 ** 3:
        return "not_hashed_large"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    tab = os.path.join(outdir, "tables")
    os.makedirs(tab, exist_ok=True)
    log = []

    def say(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        log.append(s)

    say("### A1 reporter measurement reliability")
    el, ro, paired, grp, g2b, o, qc = load()
    elements = sorted(paired.element_id)
    say(f"paired_snv elements: {len(elements)}")
    say(f"replicate_outcomes rows total: {len(ro)}; restricted to paired_snv: {len(ro[ro.element_id.isin(set(elements))])}")

    # ---- below_qc audit (full file and paired_snv restriction)
    qc_all = ro[ro.assay_state != "observed"]
    say(f"below_qc rows in full file: {len(qc_all)}; inside paired_snv elements: {len(qc)}")
    qc_tab = (qc_all.groupby(["assay_state", "missing_reason", "context_id", "allele"])
              .size().rename("rows_all_elements").reset_index())
    qc_in = (qc.groupby(["assay_state", "missing_reason", "context_id", "allele"])
             .size().rename("rows_paired_snv_elements").reset_index())
    qc_tab = qc_tab.merge(qc_in, how="left").fillna({"rows_paired_snv_elements": 0})
    qc_tab.to_csv(f"{tab}/below_qc_audit.tsv", sep="\t", index=False)
    say(qc_tab.to_string(index=False))
    qc_elems = sorted(set(qc_all.element_id))
    say(f"below_qc distinct elements (full file): {len(qc_elems)}; "
        f"distinct element-allele: {qc_all[['element_id','allele']].drop_duplicates().shape[0]}; "
        f"of those elements, paired_snv: {len(set(qc_elems) & set(elements))}")
    pd.DataFrame({"element_id": qc_elems,
                  "is_paired_snv": [e in set(elements) for e in qc_elems]}
                 ).to_csv(f"{tab}/below_qc_elements.tsv", sep="\t", index=False)

    # ---- blocks
    grp_e = grp.reindex(elements)
    blk = grp_e.map(g2b)
    n_missing_block = int(blk.isna().sum())
    say(f"elements without a long-range block: {n_missing_block}")
    blocks = sorted(blk.dropna().unique())
    say(f"distinct long-range blocks covering paired_snv elements: {len(blocks)}")
    say(f"distinct outer locus groups: {grp_e.nunique()}")
    bcode_all = pd.Series(pd.Categorical(blk, categories=blocks).codes, index=elements)

    # ---- matrices
    M = build_matrices(o, elements)

    # completeness: an element-context is usable only if all 4 replicate deltas exist
    miss_rows = []
    D = {}
    A = {}
    for ctx in CONTEXTS:
        aref = M[(ctx, "ref", "act")].values
        aalt = M[(ctx, "alt", "act")].values
        D[ctx] = aalt - aref
        A[ctx] = (aalt + aref) / 2.0
        nmiss = int(np.isnan(D[ctx]).any(axis=1).sum())
        ncell = int(np.isnan(D[ctx]).sum())
        miss_rows.append(dict(context=ctx, elements=len(elements),
                              elements_missing_at_least_one_replicate_delta=nmiss,
                              missing_element_replicate_cells=ncell))
        say(f"{ctx}: elements with >=1 missing replicate delta = {nmiss} "
            f"({ncell} missing element x replicate cells)")
    pd.DataFrame(miss_rows).to_csv(f"{tab}/missing_replicates.tsv", sep="\t", index=False)

    # ---- depth and barcode strata keys
    dna_min = {}
    bc_min = {}
    for ctx in CONTEXTS:
        dna = np.minimum(np.nanmin(M[(ctx, "ref", "DNA")].values, axis=1),
                         np.nanmin(M[(ctx, "alt", "DNA")].values, axis=1))
        bc = np.minimum(np.nanmin(M[(ctx, "ref", "n_barcodes")].values, axis=1),
                        np.nanmin(M[(ctx, "alt", "n_barcodes")].values, axis=1))
        dna_min[ctx] = dna
        bc_min[ctx] = bc

    def tertile(v, ok):
        t = np.full(len(v), -1)
        vv = v[ok]
        q1, q2 = np.nanpercentile(vv, [100 / 3, 200 / 3])
        t[ok] = np.where(vv <= q1, 0, np.where(vv <= q2, 1, 2))
        return t, (q1, q2)

    # ---- assemble the combos
    combos = {}  # key -> (H, ok_mask, stratum label)

    def add(key, Hfull, ok, strat_v=None, strat_name=None, ctx_for_key=None):
        """Register full-population combo plus its tertile strata."""
        combos[key + ("all",)] = (Hfull, ok)
        for sname, vals in (strat_v or {}).items():
            t, qs = tertile(vals, ok)
            for k in range(3):
                m = ok & (t == k)
                combos[key + (f"{sname}_tertile{k+1}",)] = (Hfull, m)

    # activity and allele effect per context
    for ctx in CONTEXTS:
        for qty, Mat in (("activity", A[ctx]), ("allele_effect", D[ctx])):
            H = np.stack([half_means(Mat, p) for p in range(3)], axis=1)
            ok = np.isfinite(H).all(axis=(1, 2)) & ~np.isnan(bcode_all.values)
            add((ctx, qty), H, ok,
                {"min_dna": dna_min[ctx], "min_barcodes": bc_min[ctx]})

    # treatment minus control allele effect, same-index convention
    for cell, (trt, ctl) in CELL_PAIRS.items():
        J = D[trt] - D[ctl]
        H = np.stack([half_means(J, p) for p in range(3)], axis=1)
        ok = np.isfinite(H).all(axis=(1, 2)) & ~np.isnan(bcode_all.values)
        add((cell, "treatment_minus_control_allele_effect"), H, ok,
            {"min_dna": np.minimum(dna_min[trt], dna_min[ctl]),
             "min_barcodes": np.minimum(bc_min[trt], bc_min[ctl])})

    say(f"combos registered: {len(combos)}")

    # ---- point estimates
    rows = []
    boot_jobs = []
    for key, (H, ok) in combos.items():
        ctx, qty, strat = key
        Hs = H[ok]
        codes = bcode_all.values[ok]
        nb_used = len(np.unique(codes))
        st = reliability_from_halves(Hs)
        base = dict(context=ctx, quantity=qty, stratum=strat,
                    n_elements=int(ok.sum()), n_blocks=int(nb_used))
        for p in range(3):
            rows.append(dict(base, partition="|".join(
                ["".join(SPLITS[p][0]), "".join(SPLITS[p][1])]),
                pearson=float(st["pearson_by_partition"][p]),
                spearman=float(st["spearman_by_partition"][p]),
                projected_4rep=sb_project(float(st["pearson_by_partition"][p])),
                projected_4rep_spearman=sb_project(float(st["spearman_by_partition"][p]))))
        rows.append(dict(base, partition="mean_of_three",
                         pearson=st["pearson"], spearman=st["spearman"],
                         projected_4rep=st["projected_4rep"],
                         projected_4rep_spearman=st["projected_4rep_spearman"]))
        # renumber blocks contiguously for the bootstrap
        uniq, codes2 = np.unique(codes, return_inverse=True)
        boot_jobs.append((key, Hs, codes2, len(uniq)))

    say(f"bootstrap jobs: {len(boot_jobs)}, {N_BOOT} draws each, seed {BOOT_SEED}")
    ncpu = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
    t0 = time.time()
    with Pool(ncpu) as pool:
        boot_out = dict(pool.map(_boot_worker, boot_jobs, chunksize=1))
    say(f"bootstrap wall seconds: {time.time() - t0:.0f}")

    R = pd.DataFrame(rows)
    for key, ci in boot_out.items():
        m = (R.context == key[0]) & (R.quantity == key[1]) & (R.stratum == key[2]) & \
            (R.partition == "mean_of_three")
        for k, v in ci.items():
            if k.startswith("_"):
                continue
            R.loc[m, k] = v
    R["bootstrap_draws"] = N_BOOT
    R["bootstrap_seed"] = BOOT_SEED
    R["resampling_unit"] = "borzoi_long_range_block"
    R.to_csv(f"{tab}/reliability_by_context.tsv", sep="\t", index=False)
    say("\n### main table (mean of three partitions, full population)")
    main_tab = R[(R.partition == "mean_of_three") & (R.stratum == "all")]
    say(main_tab[["context", "quantity", "n_elements", "n_blocks", "pearson",
                  "pearson_lo", "pearson_hi", "spearman", "projected_4rep",
                  "projected_4rep_lo", "projected_4rep_hi"]].round(4).to_string(index=False))
    say("\n### tertile strata")
    say(R[(R.partition == "mean_of_three") & (R.stratum != "all")][
        ["context", "quantity", "stratum", "n_elements", "n_blocks", "pearson",
         "projected_4rep"]].round(4).to_string(index=False))

    # ---- nine-combination sensitivity for the treatment contrast
    nine = []
    for cell, (trt, ctl) in CELL_PAIRS.items():
        Ht = np.stack([half_means(D[trt], p) for p in range(3)], axis=1)
        Hc = np.stack([half_means(D[ctl], p) for p in range(3)], axis=1)
        okc = np.isfinite(Ht).all(axis=(1, 2)) & np.isfinite(Hc).all(axis=(1, 2)) & \
            ~np.isnan(bcode_all.values)
        for pt, pc in itertools.product(range(3), range(3)):
            hA = Ht[okc, pt, 0] - Hc[okc, pc, 0]
            hB = Ht[okc, pt, 1] - Hc[okc, pc, 1]
            nine.append(dict(cell_line=cell,
                             treatment_partition="|".join(["".join(SPLITS[pt][0]),
                                                           "".join(SPLITS[pt][1])]),
                             control_partition="|".join(["".join(SPLITS[pc][0]),
                                                         "".join(SPLITS[pc][1])]),
                             same_index=bool(pt == pc),
                             n_elements=int(okc.sum()),
                             pearson=_pearson(hA, hB), spearman=_spearman(hA, hB),
                             projected_4rep=sb_project(_pearson(hA, hB))))
    N9 = pd.DataFrame(nine)
    N9.to_csv(f"{tab}/treatment_nine_combination.tsv", sep="\t", index=False)
    say("\n### nine-combination sensitivity, treatment-minus-control allele effect")
    say(N9.round(4).to_string(index=False))
    a14 = []
    for cell in CELL_PAIRS:
        s = N9[N9.cell_line == cell]
        same = float(s[s.same_index].pearson.mean())
        nine_m = float(s.pearson.mean())
        same4 = float(s[s.same_index].projected_4rep.mean())
        nine4 = float(s.projected_4rep.mean())
        a14.append(dict(cell_line=cell, same_index_pearson=same, nine_combination_pearson=nine_m,
                        abs_difference_pearson=abs(same - nine_m),
                        same_index_projected4=same4, nine_combination_projected4=nine4,
                        abs_difference_projected4=abs(same4 - nine4)))
    A14 = pd.DataFrame(a14)
    A14.to_csv(f"{tab}/a1_4_convention_difference.tsv", sep="\t", index=False)
    say(A14.round(4).to_string(index=False))

    # ---- decomposition of j: which component carries the reliable part
    # d = log2((RNA_alt+c)/(DNA_alt+c)) - log2((RNA_ref+c)/(DNA_ref+c))
    #   = [log2(RNA_alt+c) - log2(RNA_ref+c)] - [log2(DNA_alt+c) - log2(DNA_ref+c)]
    dec_rows = []
    for cell, (trt, ctl) in CELL_PAIRS.items():
        comp = {}
        for ctx in (trt, ctl):
            comp[(ctx, "rna")] = (np.log2(M[(ctx, "alt", "RNA")].values + PC)
                                  - np.log2(M[(ctx, "ref", "RNA")].values + PC))
            comp[(ctx, "dna")] = (np.log2(M[(ctx, "alt", "DNA")].values + PC)
                                  - np.log2(M[(ctx, "ref", "DNA")].values + PC))
        variants = {
            "j_full": D[trt] - D[ctl],
            "j_rna_ratio_only": comp[(trt, "rna")] - comp[(ctl, "rna")],
            "j_minus_dna_ratio_only": -(comp[(trt, "dna")] - comp[(ctl, "dna")]),
            "d_treated": D[trt],
            "d_control": D[ctl],
            "d_mean_of_conditions": (D[trt] + D[ctl]) / 2.0,
        }
        for name, Mat in variants.items():
            H = np.stack([half_means(Mat, p) for p in range(3)], axis=1)
            ok = np.isfinite(H).all(axis=(1, 2))
            st = reliability_from_halves(H[ok])
            dec_rows.append(dict(cell_line=cell, component=name, n_elements=int(ok.sum()),
                                 pearson=st["pearson"], spearman=st["spearman"],
                                 projected_4rep=st["projected_4rep"]))
        # is j just a rescaled common allele effect? cross-half (independent) correlation
        Hj = np.stack([half_means(variants["j_full"], p) for p in range(3)], axis=1)
        Hd = np.stack([half_means(variants["d_mean_of_conditions"], p) for p in range(3)], axis=1)
        Hdna = np.stack([half_means(variants["j_minus_dna_ratio_only"], p) for p in range(3)], axis=1)
        okx = (np.isfinite(Hj).all(axis=(1, 2)) & np.isfinite(Hd).all(axis=(1, 2))
               & np.isfinite(Hdna).all(axis=(1, 2)))
        cross_jd = float(np.mean([_pearson(Hj[okx, p, 0], Hd[okx, p, 1]) for p in range(3)]
                                 + [_pearson(Hj[okx, p, 1], Hd[okx, p, 0]) for p in range(3)]))
        cross_jdna = float(np.mean([_pearson(Hj[okx, p, 0], Hdna[okx, p, 1]) for p in range(3)]
                                   + [_pearson(Hj[okx, p, 1], Hdna[okx, p, 0]) for p in range(3)]))
        full_j = np.nanmean(variants["j_full"], axis=1)
        full_dna = np.nanmean(variants["j_minus_dna_ratio_only"], axis=1)
        full_rna = np.nanmean(variants["j_rna_ratio_only"], axis=1)
        full_d = np.nanmean(variants["d_mean_of_conditions"], axis=1)
        fok = np.isfinite(full_j) & np.isfinite(full_dna) & np.isfinite(full_d)
        dec_rows.append(dict(cell_line=cell, component="_cross_half_j_vs_mean_d",
                             n_elements=int(okx.sum()), pearson=cross_jd,
                             spearman=np.nan, projected_4rep=np.nan))
        dec_rows.append(dict(cell_line=cell, component="_cross_half_j_vs_minus_dna_ratio",
                             n_elements=int(okx.sum()), pearson=cross_jdna,
                             spearman=np.nan, projected_4rep=np.nan))
        dec_rows.append(dict(cell_line=cell, component="_same_half_4rep_j_vs_minus_dna_ratio",
                             n_elements=int(fok.sum()),
                             pearson=_pearson(full_j[fok], full_dna[fok]),
                             spearman=_spearman(full_j[fok], full_dna[fok]),
                             projected_4rep=np.nan))
        dec_rows.append(dict(cell_line=cell, component="_same_half_4rep_j_vs_rna_ratio",
                             n_elements=int(fok.sum()),
                             pearson=_pearson(full_j[fok], full_rna[fok]),
                             spearman=_spearman(full_j[fok], full_rna[fok]),
                             projected_4rep=np.nan))
        dec_rows.append(dict(cell_line=cell, component="_sd_of_4rep_means",
                             n_elements=int(fok.sum()),
                             pearson=float(np.std(full_j[fok], ddof=1)),
                             spearman=float(np.std(full_dna[fok], ddof=1)),
                             projected_4rep=float(np.std(full_rna[fok], ddof=1))))
        # is the reliable part of j a per-condition scale difference, j ~ beta * dbar?
        beta = float(np.polyfit(full_d[fok], full_j[fok], 1)[0])
        Mres = variants["j_full"] - beta * variants["d_mean_of_conditions"]
        Hr = np.stack([half_means(Mres, p) for p in range(3)], axis=1)
        okr = np.isfinite(Hr).all(axis=(1, 2))
        str_ = reliability_from_halves(Hr[okr])
        dec_rows.append(dict(cell_line=cell, component="j_residual_after_scale_removal",
                             n_elements=int(okr.sum()), pearson=str_["pearson"],
                             spearman=str_["spearman"],
                             projected_4rep=str_["projected_4rep"]))
        share = float(beta ** 2 * np.var(full_d[fok], ddof=1) / np.var(full_j[fok], ddof=1))
        dec_rows.append(dict(cell_line=cell, component="_scale_beta_and_variance_share",
                             n_elements=int(fok.sum()), pearson=beta, spearman=share,
                             projected_4rep=np.nan))
        for ctx in (trt, ctl):
            v = np.nanmean(D[ctx], axis=1)
            dec_rows.append(dict(cell_line=cell, component=f"_sd_4rep_mean_d_{ctx}",
                                 n_elements=int(np.isfinite(v).sum()),
                                 pearson=float(np.nanstd(v, ddof=1)),
                                 spearman=np.nan, projected_4rep=np.nan))
    DEC2 = pd.DataFrame(dec_rows)
    DEC2.to_csv(f"{tab}/treatment_contrast_decomposition.tsv", sep="\t", index=False)
    say("\n### decomposition of the treatment contrast "
        "(rows prefixed _ are diagnostics, not reliabilities; "
        "_sd_ row holds sd(j), sd(-dna ratio contrast), sd(rna ratio contrast))")
    say(DEC2.round(4).to_string(index=False))

    # ---- replicate residual correlation matrix (batch pairing diagnostic)
    samples = [(c, r) for c in CONTEXTS for r in REPS]
    names = [f"{c}_r{r}" for c, r in samples]

    def residual_matrix(get):
        X = np.column_stack([get(c)[:, REPS.index(r)] for c, r in samples])
        ok = np.isfinite(X).all(axis=1)
        X = X[ok]
        X = X - X.mean(axis=1, keepdims=True)
        X = X - X.mean(axis=0, keepdims=True)
        C = np.corrcoef(X, rowvar=False)
        return pd.DataFrame(C, index=names, columns=names), int(ok.sum())

    diag_rows = []
    rng = np.random.default_rng(PERM_SEED)
    for label, getter in (("construct_activity", lambda c: A[c]),
                          ("allele_effect", lambda c: D[c]),
                          ("log2_dna_mean_of_alleles",
                           lambda c: np.log2(
                               (M[(c, "ref", "DNA")].values + M[(c, "alt", "DNA")].values) / 2 + PC))):
        C, n_used = residual_matrix(getter)
        C.to_csv(f"{tab}/replicate_residual_correlation_{label}.tsv", sep="\t")
        say(f"\n### residual correlation matrix, {label} (n elements = {n_used})")
        say(C.round(3).to_string())
        for cell, (trt, ctl) in CELL_PAIRS.items():
            it = [names.index(f"{trt}_r{r}") for r in REPS]
            ic = [names.index(f"{ctl}_r{r}") for r in REPS]
            block = C.values[np.ix_(it, ic)]
            same = float(np.mean(np.diag(block)))
            off = float((block.sum() - np.trace(block)) / 12.0)
            obs = same - off
            null = np.empty(N_PERM)
            for b in range(N_PERM):
                perm = rng.permutation(4)
                bb = block[:, perm]
                null[b] = np.mean(np.diag(bb)) - (bb.sum() - np.trace(bb)) / 12.0
            p = float((np.sum(np.abs(null) >= abs(obs)) + 1) / (N_PERM + 1))
            diag_rows.append(dict(quantity=label, cell_line=cell,
                                  mean_same_index_cross_condition_r=same,
                                  mean_different_index_cross_condition_r=off,
                                  difference=obs, permutation_p_two_sided=p,
                                  n_permutations=N_PERM, permutation_seed=PERM_SEED))
            say(f"  {cell}: same-index cross-condition r={same:.4f} "
                f"different-index r={off:.4f} diff={obs:+.4f} perm p={p:.4f}")
        # within-condition mean r for context (DNA structure question)
        for ctx in CONTEXTS:
            ii = [names.index(f"{ctx}_r{r}") for r in REPS]
            sub = C.values[np.ix_(ii, ii)]
            within = float((sub.sum() - np.trace(sub)) / 12.0)
            diag_rows.append(dict(quantity=label, cell_line=ctx,
                                  mean_same_index_cross_condition_r=np.nan,
                                  mean_different_index_cross_condition_r=np.nan,
                                  difference=np.nan, permutation_p_two_sided=np.nan,
                                  n_permutations=0, permutation_seed=PERM_SEED,
                                  within_context_mean_r=within))
    DG = pd.DataFrame(diag_rows)
    DG.to_csv(f"{tab}/replicate_pairing_diagnostic.tsv", sep="\t", index=False)

    # ---- LX-2 DNA depth structure, explicit
    dna_rows = []
    for ctx in CONTEXTS:
        for allele in ["ref", "alt"]:
            v = M[(ctx, allele, "DNA")].values
            for i, r in enumerate(REPS):
                col = v[:, i]
                col = col[np.isfinite(col)]
                dna_rows.append(dict(context=ctx, allele=allele, replicate=r,
                                     n=len(col), total_DNA=float(col.sum()),
                                     median_DNA=float(np.median(col)),
                                     mean_log2_DNA=float(np.mean(np.log2(col + PC)))))
    DNAT = pd.DataFrame(dna_rows)
    DNAT.to_csv(f"{tab}/dna_depth_by_sample.tsv", sep="\t", index=False)
    say("\n### per-sample DNA depth")
    say(DNAT.round(3).to_string(index=False))

    # ---- reproduction of the reference 1,033-element values
    fix = pd.read_csv(FIXTURE, sep="\t", dtype=str)
    fids = sorted(set(fix.element_id))
    say(f"\n### reference reproduction on the {len(fids)}-element common fixture "
        f"({fix.outer_locus_sequence_group_id.nunique()} locus groups)")
    pos = {e: i for i, e in enumerate(elements)}
    sel = np.array([pos[e] for e in fids if e in pos])
    say(f"fixture elements present in the paired_snv set: {len(sel)} of {len(fids)}")
    ref_expect = {"LX2_TGFb": 0.9195, "LX2_control": 0.9228,
                  "HepG2_control": 0.8861, "HepG2_PAOA": 0.8580}
    rep_rows = []
    for ctx in CONTEXTS:
        Mat = D[ctx][sel]
        # (a) complete-case convention used everywhere else in this package
        H = np.stack([half_means(Mat, p) for p in range(3)], axis=1)
        ok = np.isfinite(H).all(axis=(1, 2))
        st = reliability_from_halves(H[ok])
        # (b) the reference convention: pandas .mean(axis=1) skips NaN inside a half
        Hn = np.stack([np.column_stack([
            np.nanmean(Mat[:, [REPS.index(x) for x in SPLITS[p][0]]], axis=1),
            np.nanmean(Mat[:, [REPS.index(x) for x in SPLITS[p][1]]], axis=1)])
            for p in range(3)], axis=1)
        okn = np.isfinite(Hn).all(axis=(1, 2))
        stn = reliability_from_halves(Hn[okn])
        rep_rows.append(dict(context=ctx, convention="complete_case",
                             n_elements=int(ok.sum()),
                             raw_2v2_pearson=st["pearson"],
                             projected_4rep=st["projected_4rep"],
                             ceiling_sqrt_projected=float(np.sqrt(st["projected_4rep"])),
                             reference_projected_4rep=ref_expect[ctx],
                             absolute_difference=abs(st["projected_4rep"] - ref_expect[ctx])))
        rep_rows.append(dict(context=ctx, convention="reference_nan_skipping",
                             n_elements=int(okn.sum()),
                             raw_2v2_pearson=stn["pearson"],
                             projected_4rep=stn["projected_4rep"],
                             ceiling_sqrt_projected=float(np.sqrt(stn["projected_4rep"])),
                             reference_projected_4rep=ref_expect[ctx],
                             absolute_difference=abs(stn["projected_4rep"] - ref_expect[ctx])))
    REP = pd.DataFrame(rep_rows)
    REP.to_csv(f"{tab}/reference_reproduction_1033.tsv", sep="\t", index=False)
    say(REP.round(4).to_string(index=False))

    # ---- predictions
    allele_rows = main_tab[main_tab.quantity == "allele_effect"].set_index("context")
    act_rows = main_tab[main_tab.quantity == "activity"].set_index("context")
    trt_rows = main_tab[main_tab.quantity ==
                        "treatment_minus_control_allele_effect"].set_index("context")

    preds = []
    ok11 = all(0.70 <= allele_rows.loc[c, "pearson"] <= 0.90 for c in CONTEXTS)
    preds.append(dict(prediction="A1.1",
                      statement="Raw 2v2 Pearson reliability of the allele effect is in [0.70,0.90] in all four contexts",
                      observed="; ".join(f"{c} {allele_rows.loc[c,'pearson']:.4f}" for c in CONTEXTS),
                      verdict="met" if ok11 else "not met"))
    ok12 = all(act_rows.loc[c, "pearson"] > allele_rows.loc[c, "pearson"] for c in CONTEXTS)
    preds.append(dict(prediction="A1.2",
                      statement="Activity reliability exceeds allele-effect reliability in every context",
                      observed="; ".join(
                          f"{c} activity {act_rows.loc[c,'pearson']:.4f} vs allele {allele_rows.loc[c,'pearson']:.4f}"
                          for c in CONTEXTS),
                      verdict="met" if ok12 else "not met"))
    ok13 = all(trt_rows.loc[c, "pearson"] < 0.40 for c in CELL_PAIRS)
    preds.append(dict(prediction="A1.3",
                      statement="Treatment-minus-control allele effect has raw 2v2 reliability below 0.40 in both cell lines",
                      observed="; ".join(f"{c} {trt_rows.loc[c,'pearson']:.4f}" for c in CELL_PAIRS),
                      verdict="met" if ok13 else "not met"))
    ok14 = bool((A14.abs_difference_pearson < 0.05).all())
    preds.append(dict(prediction="A1.4",
                      statement="Nine-combination convention differs from same-index by less than 0.05",
                      observed="; ".join(
                          f"{r.cell_line} |{r.same_index_pearson:.4f}-{r.nine_combination_pearson:.4f}|={r.abs_difference_pearson:.4f}"
                          for r in A14.itertuples()),
                      verdict="met" if ok14 else "not met"))
    P = pd.DataFrame(preds)
    P.to_csv(f"{tab}/predictions.tsv", sep="\t", index=False)
    say("\n### predictions")
    say(P.to_string(index=False))

    # ---- decision rule
    dec = []
    for cell in CELL_PAIRS:
        p4 = float(trt_rows.loc[cell, "projected_4rep"])
        dec.append(dict(cell_line=cell, projected_4rep=p4,
                        projected_4rep_lo=float(trt_rows.loc[cell, "projected_4rep_lo"]),
                        projected_4rep_hi=float(trt_rows.loc[cell, "projected_4rep_hi"]),
                        raw_2v2_pearson=float(trt_rows.loc[cell, "pearson"]),
                        threshold=0.50,
                        verdict=("measurable" if p4 >= 0.50 else "not measurable")))
    DEC = pd.DataFrame(dec)
    DEC.to_csv(f"{tab}/decision_rule_section3.tsv", sep="\t", index=False)
    say("\n### section 3 decision rule (projected 4-replicate reliability of j, threshold 0.50)")
    say(DEC.round(4).to_string(index=False))

    # ---- provenance
    man = []
    for p in INPUTS:
        man.append(dict(path=p, bytes=os.path.getsize(p), sha256=sha256(p)))
    pd.DataFrame(man).to_csv(f"{outdir}/MANIFEST.tsv", sep="\t", index=False)
    with open(f"{outdir}/pip_freeze.txt", "w") as fh:
        fh.write(subprocess.run([sys.executable, "-m", "pip", "freeze"],
                                capture_output=True, text=True).stdout)
    with open(f"{outdir}/seeds.json", "w") as fh:
        json.dump({"bootstrap_seed": BOOT_SEED, "bootstrap_draws": N_BOOT,
                   "permutation_seed": PERM_SEED, "permutations": N_PERM,
                   "pseudocount": PC, "python": sys.version}, fh, indent=2)
    with open(f"{outdir}/run_log.txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    say(f"\nwrote {tab}")


if __name__ == "__main__":
    main(sys.argv[1])
