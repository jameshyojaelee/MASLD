#!/usr/bin/env python3
"""Robustness of the coding-mass asymmetry: the three checks HIGH_IMPACT row 7
records as still owed, plus an LD-block clustering arm.

The headline compares the share of fine-mapped posterior mass in protein-altering
or canonical-splice variants between direct-disease and liver-enzyme credible
sets. After the PIP-sum gate (2026-09-23 correction) it is 22.68% against 3.39%,
ratio 6.69, over 207 direct and 6,208 enzyme credible-set instances in 24 and 259
operational clusters. Those clusters are 264 physical 1-Mb regions; 19 hold both
classes. Three things could produce it without the trait class being the cause:

  1  ONE STUDY. MVP_NAFLD_AMR contributes 95 of the direct instances. If the
     result is that study, it is a study result, not a trait-class result.
  2  ONE CLUSTER. 24 clusters is few enough that a single coding-heavy locus
     could carry the fraction.
  3  THE WINDOW. 1 Mb is arbitrary; a different window changes how many clusters
     exist and therefore what is being averaged.

A fourth arm replaces the arbitrary window with the fine-mapping run's LD blocks,
so clusters are defined by measured LD rather than by distance.

Everything reuses the producer's functions, so the gate, the physical-cluster
bootstrap and the estimand are identical to `fig2_multiplicity_correction.py`;
only the subset or the cluster definition changes. As there, no permutation p is
computed. Each fit draws from its own generator seeded by SEED and its fit id,
so a fit's interval does not depend on the order the fits run in.
"""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig2_multiplicity_correction import (            # noqa: E402
    SCOPES, CODING_COLS, ARCH_COLS, SEED, N_BOOT, EXPECTED_GATE_FAILURES,
    CORRECTED_CODING, apply_gate, fit_rng, parse_locus, assign_clusters,
    cluster_table, analyse, sha256, bh)

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
    """Cluster by the fine-mapping run's own LD blocks, looked up in the table of
    each set's ancestry. The five ancestry tables share one set of block
    coordinates (asserted below), so the cluster id is the block coordinates
    alone: the same region reached from two ancestries or two trait classes is
    one physical cluster. A credible set that falls in no block keeps a
    singleton id, which is conservative: it cannot merge with anything."""
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
            blocks[anc] = b[["chr", "start", "stop"]].reset_index(drop=True)
    coords = list(blocks.values())
    if coords and not all(c.equals(coords[0]) for c in coords[1:]):
        raise ValueError("ancestry LD-block tables differ; an ancestry-free block id is invalid")
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
                cid = f"{h.chr}:{h.start}-{h.stop}"
        ids.append(cid if cid else f"singleton:{i}")
    f["cluster_id"] = ids
    return f


def run(tab, **tags):
    fit_id = f"{tags['arm']}|{tags['dropped']}"
    r = analyse(tab, fit_rng(fit_id))
    r.update(tags, fit_id=fit_id)
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    root, out = Path(args.project_root), Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    base = root / ("GWAS/finemapping/results/candidates/"
                   "noncoding-dna-precoloc-2026-08-12/run_19783321")
    arch_path = base / "credible_set_pip_architecture.tsv"
    members_path = base / "credible_set_members_annotated.tsv"
    arch = pd.read_csv(arch_path, sep="\t")
    members = pd.read_csv(members_path, sep="\t", low_memory=False)
    arch = arch[arch["trait_scope"].isin(SCOPES)].copy()
    members = members[members["trait_scope"].isin(SCOPES)].copy()
    # same credible-set gate and exclusions as the producer
    arch, _, excluded = apply_gate(arch, members)
    excluded.to_csv(out / "excluded_credible_sets.tsv", sep="\t", index=False)
    gate_fail = (excluded[excluded["reason"] == "pip_sum_gate_failed"]
                 .groupby("trait_scope").size().to_dict())
    assert gate_fail == EXPECTED_GATE_FAILURES, gate_fail
    arch = coding_frame(parse_locus(arch))

    def ci(r):
        return (f"ratio {r.get('ratio', float('nan')):.3f} "
                f"[{r.get('ratio_boot_lo', float('nan')):.3f}, {r.get('ratio_boot_hi', float('nan')):.3f}] "
                f"diff lo {r.get('diff_boot_lo', float('nan')):.4f}")

    # ---- reference: reproduce the corrected 1-Mb result before varying anything
    ref = assign_clusters(arch.copy())
    r0 = run(cluster_table(ref, "coding_mass", "arch_total"),
             arm="reference_1Mb", dropped="none")
    print(f"[ref] direct {r0['fraction_direct']:.4f} enzyme {r0['fraction_enzyme']:.4f} "
          f"{ci(r0)} physical {r0['n_physical_clusters']} shared {r0['n_shared_clusters']}",
          flush=True)
    for key, expected in CORRECTED_CODING.items():
        tol = 1e-9 if key.startswith("fraction") else 1e-7 * expected
        assert abs(r0[key] - expected) <= tol, (key, r0[key], expected)
    assert (r0["n_physical_clusters"], r0["n_shared_clusters"]) == (264, 19), r0
    print("[ref] corrected gated estimate reproduced", flush=True)

    rows = [r0]

    # ---- 1. leave-one-study-out, both classes
    for scope, tag in ((SCOPES[0], "direct"), (SCOPES[1], "enzyme")):
        for st in sorted(arch.loc[arch.trait_scope == scope, "study"].unique()):
            sub = assign_clusters(arch[~((arch.trait_scope == scope) & (arch.study == st))].copy())
            rows.append(run(cluster_table(sub, "coding_mass", "arch_total"),
                            arm=f"leave_one_study_out_{tag}", dropped=st,
                            n_instances_dropped=int(((arch.trait_scope == scope) &
                                                     (arch.study == st)).sum())))
            print(f"[loso-{tag}] -{st:<32} {ci(rows[-1])}", flush=True)

    # ---- 2. leave-one-cluster-out, direct clusters only (24 of them)
    for cid in sorted(ref.loc[ref.trait_scope == SCOPES[0], "cluster_id"].unique()):
        sub = ref[~((ref.trait_scope == SCOPES[0]) & (ref.cluster_id == cid))].copy()
        rows.append(run(cluster_table(sub, "coding_mass", "arch_total"),
                        arm="leave_one_cluster_out_direct", dropped=cid))
    print(f"[loco] {sum(r['arm']=='leave_one_cluster_out_direct' for r in rows)} direct "
          f"clusters dropped one at a time", flush=True)

    # ---- 3. window size
    for gap in (250_000, 500_000, 1_000_000, 2_000_000, 5_000_000):
        sub = assign_clusters_gap(arch.copy(), gap)
        rows.append(run(cluster_table(sub, "coding_mass", "arch_total"),
                        arm="window_size", dropped=f"{gap/1e6:g}Mb"))
        print(f"[window] {gap/1e6:g} Mb  direct clusters "
              f"{rows[-1]['n_clusters_direct']} enzyme {rows[-1]['n_clusters_enzyme']} "
              f"physical {rows[-1].get('n_physical_clusters')} "
              f"shared {rows[-1].get('n_shared_clusters')} {ci(rows[-1])}", flush=True)

    # ---- 4. LD-block clustering, one physical block id across ancestries
    ld = assign_clusters_ldblock(arch.copy(), root)
    rows.append(run(cluster_table(ld, "coding_mass", "arch_total"),
                    arm="ld_block_clustering", dropped="none"))
    nsing = int(ld.cluster_id.str.startswith("singleton:").sum())
    rows[-1]["n_singletons_no_block"] = nsing
    print(f"[ldblock] direct clusters {rows[-1]['n_clusters_direct']} "
          f"enzyme {rows[-1]['n_clusters_enzyme']} "
          f"physical {rows[-1].get('n_physical_clusters')} "
          f"shared {rows[-1].get('n_shared_clusters')} {ci(rows[-1])} "
          f"({nsing} sets fell in no block)", flush=True)

    assert len(rows) == 62, len(rows)
    res = pd.DataFrame(rows)
    res["seed"] = SEED
    res["n_boot"] = N_BOOT
    res.to_csv(out / "coding_mass_robustness.tsv", sep="\t", index=False)

    # ---- summary the manuscript can quote
    def rng_of(arm):
        s = res[res.arm == arm]["ratio"].dropna()
        return (float(s.min()), float(s.max()), int(len(s))) if len(s) else (np.nan, np.nan, 0)
    summ = {"inference_unit": "physical_cluster_pooled_over_trait_classes",
            "p_permutation": "not computed; see fig2_multiplicity_correction.py docstring",
            "reference_ratio_1Mb": float(r0["ratio"]),
            "reference_ratio_boot_ci": [float(r0["ratio_boot_lo"]), float(r0["ratio_boot_hi"])],
            "reference_diff_boot_ci": [float(r0["diff_boot_lo"]), float(r0["diff_boot_hi"])],
            # percentile-inversion diagnostic floored at 2/(n+1), not a test
            "reference_nominal_p_bootstrap_diff": float(r0["nominal_p_bootstrap_diff"])}
    for arm in ("leave_one_study_out_direct", "leave_one_study_out_enzyme",
                "leave_one_cluster_out_direct", "window_size"):
        lo, hi, n = rng_of(arm)
        sub = res[res.arm == arm]
        summ[arm] = {"n_fits": n, "ratio_min": lo, "ratio_max": hi,
                     "min_ratio_boot_lo": float(sub["ratio_boot_lo"].min()) if n else None,
                     "min_diff_boot_lo": float(sub["diff_boot_lo"].min()) if n else None,
                     "max_nominal_p_bootstrap_diff":
                         float(sub["nominal_p_bootstrap_diff"].max()) if n else None}
    ldr = res[res.arm == "ld_block_clustering"].iloc[0]
    summ["ld_block_clustering"] = {"ratio": float(ldr["ratio"]),
                                   "ratio_boot_ci": [float(ldr["ratio_boot_lo"]),
                                                     float(ldr["ratio_boot_hi"])],
                                   "diff_boot_ci": [float(ldr["diff_boot_lo"]),
                                                    float(ldr["diff_boot_hi"])],
                                   "nominal_p_bootstrap_diff":
                                       float(ldr["nominal_p_bootstrap_diff"]),
                                   "n_clusters_direct": int(ldr["n_clusters_direct"]),
                                   "n_clusters_enzyme": int(ldr["n_clusters_enzyme"]),
                                   "n_physical_clusters": int(ldr["n_physical_clusters"]),
                                   "n_shared_clusters": int(ldr["n_shared_clusters"])}
    most = res[res.arm == "leave_one_study_out_direct"].nsmallest(1, "ratio")
    if len(most):
        summ["most_influential_direct_study"] = {
            "study": str(most.iloc[0]["dropped"]),
            "ratio_without_it": float(most.iloc[0]["ratio"]),
            "ratio_boot_ci_without_it": [float(most.iloc[0]["ratio_boot_lo"]),
                                         float(most.iloc[0]["ratio_boot_hi"])],
            "diff_boot_ci_without_it": [float(most.iloc[0]["diff_boot_lo"]),
                                        float(most.iloc[0]["diff_boot_hi"])]}
    json.dump(summ, open(out / "robustness_summary.json", "w"), indent=1)

    (out / "input_manifest.json").write_text(json.dumps(
        [{"role": r, "path": str(p), "sha256": sha256(p)}
         for r, p in [("architecture", arch_path), ("members", members_path)]], indent=2))
    (out / "environment.txt").write_text(
        f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
        f"seed {SEED} (per fit: fit_rng(arm|dropped))\nn_boot {N_BOOT}\n"
        + subprocess.run([sys.executable, "-m", "pip", "freeze"],
                         capture_output=True, text=True).stdout)
    print("\n=== SUMMARY ===")
    print(json.dumps(summ, indent=1))
    print(f"[robust] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
