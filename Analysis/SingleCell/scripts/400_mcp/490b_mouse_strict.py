#!/usr/bin/env python
"""
490b_mouse_strict.py — Strict 1:1 ortholog-only mouse projection (no case-conversion fallback).
Addresses reviewer concern that case-conversion can false-map paralogs.

--null mode (Echo Task 19 / A15): for each program, draw N (=program size in mouse-strict
space) genes from the global atlas HVG pool 1000 times, score each random gene set on the
mouse atlas via sc.tl.score_genes, and compare the observed program mean mouse score to
the null distribution. Reports per-program permutation p / q.
"""
from __future__ import annotations
import argparse, os, gzip, re
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
ATLAS = ROOT / "Analysis/Deconvolution/reference/reference_mouse.h5ad"
ORTHO = ROOT / "archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
GENCODE = ROOT / "data/gencode_v49_gene_metadata.tsv.gz"
HVG_FILE = MCP / "cnmf_runs/global/global.overdispersed_genes.txt"
OUT = MCP / "validation/mouse"


def _build_chain_maps():
    gm = pd.read_csv(GENCODE, sep="\t")
    sym2ens = dict(zip(gm["gene_name"], gm["ensembl_base"]))
    o = pd.read_csv(ORTHO, sep="\t")
    o = o[o["orthology_type"] == "ortholog_one2one"]
    h2m = dict(zip(o["human_ensembl_gene_id"], o["mouse_ensembl_gene_id"]))

    gtf = ROOT / "../../../../home/jameslee/reference_genome/refdata-gex-GRCm39-2024-A/annotation/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz"
    mouse_gtf = Path(str(gtf))
    m_ens2sym = {}
    if mouse_gtf.exists():
        pat = re.compile(r'gene_id "([^"]+)".*gene_name "([^"]+)"')
        with gzip.open(mouse_gtf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if "\tgene\t" not in line:
                    continue
                m = pat.search(line)
                if m:
                    m_ens2sym[m.group(1).split(".")[0]] = m.group(2)
    return sym2ens, h2m, m_ens2sym


def _strict_map(human_syms, sym2ens, h2m, m_ens2sym, mouse_var_set):
    out = []
    for hs in human_syms:
        e = sym2ens.get(hs)
        if not e:
            continue
        mens = h2m.get(e)
        if not mens:
            continue
        msym = m_ens2sym.get(mens)
        if not msym:
            continue
        if msym in mouse_var_set:
            out.append(msym)
    return out


def main(k: int = 16):
    topg = pd.read_csv(MCP / f"cnmf_annot/global/program_topgenes.k{k}.tsv", sep="\t")
    sym2ens, h2m, m_ens2sym = _build_chain_maps()
    print(f"[490b] ensembl→symbol map size: {len(m_ens2sym)}")

    a = sc.read_h5ad(ATLAS)
    a.obs_names_make_unique()
    if "log1p" not in a.uns:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    mouse_var_set = set(a.var_names)

    rows = []
    for p, g in topg.groupby("program"):
        h_syms = g["gene_name"].tolist()[:100]
        strict = _strict_map(h_syms, sym2ens, h2m, m_ens2sym, mouse_var_set)
        rows.append({"program": str(p), "n_human_top": len(h_syms), "n_strict_1to1": len(strict)})
        if len(strict) >= 5:
            sc.tl.score_genes(a, strict, score_name=f"prog_{p}_strict")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"mouse_projection_strict_k{k}.tsv", sep="\t", index=False)
    print(f"[490b] Strict 1:1 mapping: mean {df['n_strict_1to1'].mean():.1f}/100 genes per program")
    print(df.to_string())

    prog_cols = [c for c in a.obs.columns if c.endswith("_strict")]
    if prog_cols:
        agg = a.obs.groupby("sample", observed=True)[prog_cols].mean()
        agg.to_csv(OUT / f"mouse_projection_strict_k{k}_mean_by_sample.tsv", sep="\t")


def null_mode(k: int = 16, n_perm: int = 1000, seed: int = 42):
    """Permutation null: random HVG draws → strict ortholog chain → mouse score.

    Per program, draw `match_to_size` random human HVG every iteration (size = number of
    one-to-one mouse-mappable genes for that program), pass through the same strict chain,
    score on the mouse atlas. Compare observed program mean mouse score to null.
    """
    rng = np.random.default_rng(seed)

    topg = pd.read_csv(MCP / f"cnmf_annot/global/program_topgenes.k{k}.tsv", sep="\t")
    sym2ens, h2m, m_ens2sym = _build_chain_maps()

    if not HVG_FILE.exists():
        raise FileNotFoundError(f"HVG file missing: {HVG_FILE}")
    hvg_pool = [line.strip() for line in HVG_FILE.read_text().splitlines() if line.strip()]
    print(f"[490b.null] global HVG pool size: {len(hvg_pool)}")

    a = sc.read_h5ad(ATLAS)
    a.obs_names_make_unique()
    if "log1p" not in a.uns:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    mouse_var_set = set(a.var_names)

    # Score each program (observed)
    obs_scores = {}
    prog_sizes = {}
    for p, g in topg.groupby("program"):
        h_syms = g["gene_name"].tolist()[:100]
        strict = _strict_map(h_syms, sym2ens, h2m, m_ens2sym, mouse_var_set)
        prog_sizes[str(p)] = len(strict)
        if len(strict) < 5:
            obs_scores[str(p)] = np.nan
            continue
        sc.tl.score_genes(a, strict, score_name=f"obs_{p}_strict")
        obs_scores[str(p)] = float(a.obs[f"obs_{p}_strict"].mean())

    # Build null distributions program-by-program: each program needs its own size match.
    # To save compute, batch by program-size: many programs share similar sizes.
    sizes_unique = sorted({s for s in prog_sizes.values() if s >= 5})
    print(f"[490b.null] unique program sizes: {sizes_unique}")

    null_means_by_size: dict[int, np.ndarray] = {}
    for sz in sizes_unique:
        means = np.full(n_perm, np.nan, dtype=float)
        for it in range(n_perm):
            # draw HVG, push through chain, target sz mouse-mappable genes
            tries = 0
            while True:
                tries += 1
                draw = rng.choice(hvg_pool, size=min(len(hvg_pool), sz * 6 + 50), replace=False)
                strict = _strict_map(list(draw), sym2ens, h2m, m_ens2sym, mouse_var_set)
                if len(strict) >= sz:
                    strict = strict[:sz]
                    break
                if tries > 5:
                    # fallback: take whatever we got (avoids infinite loop)
                    if len(strict) >= 5:
                        break
                    strict = []
                    break
            if not strict or len(strict) < 5:
                means[it] = np.nan
                continue
            tmpcol = f"null_{sz}_{it}"
            try:
                sc.tl.score_genes(a, strict, score_name=tmpcol)
                means[it] = float(a.obs[tmpcol].mean())
            except Exception:
                means[it] = np.nan
            finally:
                if tmpcol in a.obs.columns:
                    a.obs.drop(columns=[tmpcol], inplace=True)
            if (it + 1) % 100 == 0:
                print(f"[490b.null]   size={sz} iter={it + 1}/{n_perm}")
        null_means_by_size[sz] = means
        print(f"[490b.null] size={sz}: null mean={np.nanmean(means):.4f} sd={np.nanstd(means):.4f}")

    # Per-program p / q against its size-matched null
    rows = []
    for p, sz in prog_sizes.items():
        obs = obs_scores[p]
        if sz < 5 or np.isnan(obs):
            rows.append({
                "program": p, "n_strict": sz,
                "obs_mean_score": obs,
                "null_mean": np.nan, "null_sd": np.nan,
                "perm_p_two_sided": np.nan,
                "z_score": np.nan,
                "n_perm_finite": 0,
            })
            continue
        nulls = null_means_by_size[sz]
        nulls = nulls[np.isfinite(nulls)]
        n_finite = len(nulls)
        if n_finite < 10:
            rows.append({
                "program": p, "n_strict": sz, "obs_mean_score": obs,
                "null_mean": np.nan, "null_sd": np.nan,
                "perm_p_two_sided": np.nan, "z_score": np.nan,
                "n_perm_finite": n_finite,
            })
            continue
        nm = float(nulls.mean())
        ns = float(nulls.std(ddof=1)) if n_finite > 1 else np.nan
        # Two-sided permutation p with +1 smoothing
        more_extreme = (np.abs(nulls - nm) >= abs(obs - nm)).sum()
        p_perm = (more_extreme + 1) / (n_finite + 1)
        z = (obs - nm) / ns if (ns and ns > 0) else np.nan
        rows.append({
            "program": p, "n_strict": sz, "obs_mean_score": obs,
            "null_mean": nm, "null_sd": ns,
            "perm_p_two_sided": p_perm, "z_score": z,
            "n_perm_finite": n_finite,
        })

    df = pd.DataFrame(rows).sort_values("program")
    # BH FDR over programs with finite p
    pmask = df["perm_p_two_sided"].notna()
    if pmask.any():
        from scipy.stats import false_discovery_control
        try:
            df.loc[pmask, "perm_q_bh"] = false_discovery_control(
                df.loc[pmask, "perm_p_two_sided"].values, method="bh"
            )
        except Exception:
            # Fallback manual BH
            pv = df.loc[pmask, "perm_p_two_sided"].values
            order = np.argsort(pv)
            ranks = np.empty_like(order)
            ranks[order] = np.arange(1, len(pv) + 1)
            qv = pv * len(pv) / ranks
            qv = np.minimum.accumulate(qv[order[::-1]])[::-1]
            qfull = np.empty_like(pv)
            qfull[order] = qv
            df.loc[pmask, "perm_q_bh"] = qfull

    OUT.mkdir(parents=True, exist_ok=True)
    out_f = OUT / f"mouse_conservation_null_k{k}.tsv"
    df.to_csv(out_f, sep="\t", index=False)
    print(f"[490b.null] wrote {out_f}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--null", action="store_true", help="Run permutation null mode (Echo A15)")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    if args.null:
        null_mode(k=args.k, n_perm=args.n_perm, seed=args.seed)
    else:
        main(args.k)
