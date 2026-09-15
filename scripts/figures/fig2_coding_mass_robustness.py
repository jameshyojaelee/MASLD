#!/usr/bin/env python3
"""Robustness of the coding-mass asymmetry: the three checks HIGH_IMPACT row 7
records as still owed, plus an LD-block clustering arm.

The headline (direct-disease credible sets carry 22.5% of their fine-mapped
posterior mass in protein-altering or canonical-splice variants against 3.4% for
liver-enzyme sets, ratio 6.55, permutation BH q = 6.0e-4) rests on a 1-Mb
single-linkage clustering of 208 direct and 6,259 enzyme credible-set instances
into 24 and 266 operational clusters. Three things could produce it without the
trait class being the cause:

  1  ONE STUDY. MVP_NAFLD_AMR contributes 95 of the 208 direct instances. If the
     result is that study, it is a study result, not a trait-class result.
  2  ONE CLUSTER. 24 clusters is few enough that a single coding-heavy locus
     could carry the fraction.
  3  THE WINDOW. 1 Mb is arbitrary; a different window changes how many clusters
     exist and therefore what is being averaged.

A fourth arm replaces the arbitrary window with the fine-mapping run's own
ancestry-matched LD blocks, so clusters are defined by measured LD rather than
by distance.

Everything reuses the adopted producer's functions, so the clustering, the
estimand and the permutation are identical to the published test; only the
subset or the cluster definition changes.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig2_multiplicity_correction import (            # noqa: E402
    SCOPES, CODING_COLS, ARCH_COLS, SEED, N_PERM,
    parse_locus, assign_clusters, cluster_table, analyse, sha256, bh)

BLOCKS = "GWAS/finemapping/runs/uniform35_v3_2026-08-03/config/ld_blocks"


def coding_frame(arch: pd.DataFrame) -> pd.DataFrame:
    a = arch.copy()
    a["coding_mass"] = a[CODING_COLS].sum(axis=1)
    a["arch_total"] = a[ARCH_COLS].sum(axis=1)
    return a


def assign_clusters_gap(frame: pd.DataFrame, gap_bp: int) -> pd.DataFrame:
    """assign_clusters with an explicit gap. Mirrors the producer exactly."""
    f = frame.sort_values(["chrom","pos"]).copy()
    keys = []
    for chrom, block in f.groupby("chrom", sort=False):
        block = block.sort_values("pos")
        pos = block["pos"].to_numpy()
        new = np.concatenate([[True], np.diff(pos) > gap_bp])
        keys.append(pd.Series([f"{chrom}:{i}" for i in np.cumsum(new)], index=block.index))
    f["cluster_id"] = pd.concat(keys).reindex(f.index)
    return f


def assign_clusters_ldblock(frame: pd.DataFrame, root: Path) -> pd.DataFrame:
    """Cluster by the fine-mapping run's own LD blocks, ancestry-matched.
    A credible set that falls in no block for its ancestry keeps a singleton id,
    which is conservative: it cannot merge with anything."""
    blocks = {}
    for anc in ("AFR", "AMR", "EAS", "EUR", "SAS"):
        fp = root / BLOCKS / f"{anc}.tsv"
        if fp.exists():
            b = pd.read_csv(fp, sep="\t")
            cols = {c.lower(): c for c in b.columns}
            b = b.rename(columns={cols.get("chr", cols.get("chromosome", "chr")): "chr",
                                  cols.get("start", cols.get("block_start", "start")): "start",
                                  cols.get("stop", cols.get("block_stop", "stop")): "stop"})
            b["chr"] = b["chr"].astype(str).str.replace("chr", "", regex=False)
            blocks[anc] = b[["chr", "start", "stop"]]
    f = frame.copy()
    ids = []
    for i, r in f.iterrows():
        anc = r.get("ancestry", "EUR")
        b = blocks.get(anc)
        cid = None
        if b is not None:
            hit = b[(b.chr == str(r["chrom"])) &
                    (b.start <= r["pos"]) & (b.stop >= r["pos"])]
            if len(hit):
                h = hit.iloc[0]
                cid = f"{anc}:{h.chr}:{h.start}-{h.stop}"
        ids.append(cid if cid else f"singleton:{i}")
    f["cluster_id"] = ids
    return f


def run(tab, rng, **tags):
    r = analyse(tab, rng)
    r.update(tags)
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    root, out = Path(args.project_root), Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    base = root / ("GWAS/finemapping/results/candidates/"
                   "noncoding-dna-precoloc-2026-08-12/run_19783321")
    arch_path = base / "credible_set_pip_architecture.tsv"
    arch = pd.read_csv(arch_path, sep="\t")
    arch = arch[arch["trait_scope"].isin(SCOPES)].copy()
    arch = coding_frame(parse_locus(arch))

    # ---- reference: reproduce the adopted 1-Mb result before varying anything
    ref = assign_clusters(arch.copy())
    r0 = run(cluster_table(ref, "coding_mass", "arch_total"), rng,
             arm="reference_1Mb", dropped="none")
    print(f"[ref] direct {r0['fraction_direct']:.4f} enzyme {r0['fraction_enzyme']:.4f} "
          f"ratio {r0['ratio']:.3f} p_perm {r0['p_permutation']:.4g}", flush=True)
    assert abs(r0["ratio"] - 6.549352416207798) < 1e-6, r0["ratio"]
    print("[ref] adopted ratio 6.5494 reproduced exactly", flush=True)

    rows = [r0]

    # ---- 1. leave-one-study-out, both classes
    for scope, tag in ((SCOPES[0], "direct"), (SCOPES[1], "enzyme")):
        for st in sorted(arch.loc[arch.trait_scope == scope, "study"].unique()):
            sub = assign_clusters(arch[~((arch.trait_scope == scope) & (arch.study == st))].copy())
            rows.append(run(cluster_table(sub, "coding_mass", "arch_total"), rng,
                            arm=f"leave_one_study_out_{tag}", dropped=st,
                            n_instances_dropped=int(((arch.trait_scope == scope) &
                                                     (arch.study == st)).sum())))
            print(f"[loso-{tag}] -{st:<32} ratio {rows[-1].get('ratio', float('nan')):.3f} "
                  f"p {rows[-1].get('p_permutation', float('nan')):.4g}", flush=True)

    # ---- 2. leave-one-cluster-out, direct clusters only (24 of them)
    for cid in sorted(ref.loc[ref.trait_scope == SCOPES[0], "cluster_id"].unique()):
        sub = ref[~((ref.trait_scope == SCOPES[0]) & (ref.cluster_id == cid))].copy()
        rows.append(run(cluster_table(sub, "coding_mass", "arch_total"), rng,
                        arm="leave_one_cluster_out_direct", dropped=cid))
    print(f"[loco] {sum(r['arm']=='leave_one_cluster_out_direct' for r in rows)} direct "
          f"clusters dropped one at a time", flush=True)

    # ---- 3. window size
    for gap in (250_000, 500_000, 1_000_000, 2_000_000, 5_000_000):
        sub = assign_clusters_gap(arch.copy(), gap)
        rows.append(run(cluster_table(sub, "coding_mass", "arch_total"), rng,
                        arm="window_size", dropped=f"{gap/1e6:g}Mb"))
        print(f"[window] {gap/1e6:g} Mb  direct clusters "
              f"{rows[-1]['n_clusters_direct']} enzyme {rows[-1]['n_clusters_enzyme']} "
              f"ratio {rows[-1].get('ratio', float('nan')):.3f} "
              f"p {rows[-1].get('p_permutation', float('nan')):.4g}", flush=True)

    # ---- 4. LD-block clustering, ancestry-matched
    ld = assign_clusters_ldblock(arch.copy(), root)
    rows.append(run(cluster_table(ld, "coding_mass", "arch_total"), rng,
                    arm="ld_block_clustering", dropped="none"))
    nsing = int(ld.cluster_id.str.startswith("singleton:").sum())
    rows[-1]["n_singletons_no_block"] = nsing
    print(f"[ldblock] direct clusters {rows[-1]['n_clusters_direct']} "
          f"enzyme {rows[-1]['n_clusters_enzyme']} ratio {rows[-1].get('ratio', float('nan')):.3f} "
          f"p {rows[-1].get('p_permutation', float('nan')):.4g} "
          f"({nsing} sets fell in no ancestry-matched block)", flush=True)

    res = pd.DataFrame(rows)
    res["seed"] = SEED
    res["n_perm"] = N_PERM
    res.to_csv(out / "coding_mass_robustness.tsv", sep="\t", index=False)

    # ---- summary the manuscript can quote
    def rng_of(arm):
        s = res[res.arm == arm]["ratio"].dropna()
        return (float(s.min()), float(s.max()), int(len(s))) if len(s) else (np.nan, np.nan, 0)
    summ = {"reference_ratio_1Mb": float(r0["ratio"]),
            "reference_p_permutation": float(r0["p_permutation"])}
    for arm in ("leave_one_study_out_direct", "leave_one_study_out_enzyme",
                "leave_one_cluster_out_direct", "window_size"):
        lo, hi, n = rng_of(arm)
        summ[arm] = {"n_fits": n, "ratio_min": lo, "ratio_max": hi,
                     "max_p_permutation": float(res[res.arm == arm]["p_permutation"].max())
                                          if n else None,
                     "all_p_below_0.05": bool((res[res.arm == arm]["p_permutation"] < 0.05).all())
                                          if n else None}
    ldr = res[res.arm == "ld_block_clustering"].iloc[0]
    summ["ld_block_clustering"] = {"ratio": float(ldr["ratio"]),
                                   "p_permutation": float(ldr["p_permutation"]),
                                   "n_clusters_direct": int(ldr["n_clusters_direct"]),
                                   "n_clusters_enzyme": int(ldr["n_clusters_enzyme"])}
    most = res[res.arm == "leave_one_study_out_direct"].nsmallest(1, "ratio")
    if len(most):
        summ["most_influential_direct_study"] = {
            "study": str(most.iloc[0]["dropped"]),
            "ratio_without_it": float(most.iloc[0]["ratio"]),
            "p_without_it": float(most.iloc[0]["p_permutation"])}
    json.dump(summ, open(out / "robustness_summary.json", "w"), indent=1)

    (out / "input_manifest.json").write_text(json.dumps(
        [{"role": "architecture", "path": str(arch_path), "sha256": sha256(arch_path)}], indent=2))
    (out / "environment.txt").write_text(
        f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
        f"seed {SEED}\nn_perm {N_PERM}\n")
    print("\n=== SUMMARY ===")
    print(json.dumps(summ, indent=1))
    print(f"[robust] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
