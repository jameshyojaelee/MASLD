#!/usr/bin/env python
"""Family-wise multiplicity correction for the Figure 2 trait-directness tests.

WHY. Six direct-vs-enzyme comparisons are made on the same set of operational
locus clusters: one on coding PIP mass and five on noncoding regulatory context
(promoter-proximal, ABC enhancer, lineage-accessible, resolved-no-context,
unresolved). Each is reported with its own uncorrected p. No correction is
applied across them. The family is the six tests, and it must be named and
corrected before any of them is called significant.

This matters because the coding claim's permutation p is 0.055
(`ascertainment_control/20260817T152705Z/`), so it is already at the boundary
before any correction.

TWO NULLS, BOTH REPORTED. `fig2_trait_directness_stats.py` reports
`p_two_sided_diff_gt0 = 2 * min((diff <= 0).mean(), (diff >= 0).mean())`, a
BOOTSTRAP p asking whether the direct-minus-enzyme difference is distinguishable
from zero given within-arm sampling variability. This script reproduces that and
adds a PERMUTATION p that shuffles the trait-class label across locus clusters,
asking whether the grouping itself carries the difference. The permutation null
is the appropriate test for a group-contrast claim and is the one BH is applied
to here. The bootstrap p is carried alongside so the change is auditable.

UNIT AND ITS NAME. The unit is a 1-Mb single-linkage cluster of credible-set
sentinel positions. These are OPERATIONAL CLUSTERS, not independently segregating
loci: no LD was consulted in forming them, a 1-Mb window can span more than one
independent signal, and long-range LD can link clusters. Every output here uses
`operational cluster`. The estimand within a cluster is unchanged from the
published script (sum of per-cluster mean masses over sum of per-cluster mean
denominators).

Read-only inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260817
N_BOOT = 10000
N_PERM = 10000
LOCUS_GAP_BP = 1_000_000

SCOPES = ["tier1_direct_masld_pdff", "tier2_liver_enzyme"]
CODING_COLS = ["protein_altering_pip_mass", "canonical_splice_pip_mass"]
ARCH_COLS = CODING_COLS + ["synonymous_or_utr_pip_mass", "other_noncoding_pip_mass"]
CTX_COLS = ["promoter_mass", "abc_mass", "accessible_mass",
            "resolved_no_context_mass", "unresolved_context_mass"]

TEST_LABEL = {
    "coding_mass": "Coding PIP mass (protein-altering + canonical splice)",
    "promoter_mass": "Promoter-proximal (noncoding)",
    "abc_mass": "ABC enhancer overlap (noncoding)",
    "accessible_mass": "Lineage-accessible (noncoding)",
    "resolved_no_context_mass": "Resolved, no context (noncoding)",
    "unresolved_context_mass": "Unresolved context (noncoding)",
}

PUBLISHED = {
    ("tier1_direct_masld_pdff", "protein_altering_pip_mass"): 0.157736519057125,
    ("tier2_liver_enzyme", "protein_altering_pip_mass"): 0.0325683421565398,
    ("tier1_direct_masld_pdff", "canonical_splice_pip_mass"): 0.000361875951271015,
    ("tier2_liver_enzyme", "canonical_splice_pip_mass"): 0.0,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def parse_locus(frame: pd.DataFrame) -> pd.DataFrame:
    """Copied verbatim from fig2_trait_directness_stats.py."""
    text = frame["locus"].astype(str)
    frame = frame.copy()
    frame["chrom"] = text.str.split(".").str[0]
    frame["pos"] = pd.to_numeric(text.str.split(".").str[1], errors="coerce")
    if frame["pos"].isna().any():
        raise ValueError("unparsable locus positions present")
    return frame


def assign_clusters(frame: pd.DataFrame) -> pd.DataFrame:
    """Copied verbatim from fig2_trait_directness_stats.py (`assign_independent_loci`),
    renamed because the output is an operational cluster, not an independent locus."""
    frame = frame.sort_values(["chrom", "pos"]).copy()
    keys = []
    for chrom, block in frame.groupby("chrom", sort=False):
        pos = block["pos"].to_numpy()
        new_cluster = np.concatenate([[True], np.diff(pos) > LOCUS_GAP_BP])
        keys.append(pd.Series(
            [f"{chrom}:{i}" for i in np.cumsum(new_cluster)], index=block.index))
    frame["cluster_id"] = pd.concat(keys).reindex(frame.index)
    return frame


def build_context_masses(members: pd.DataFrame) -> pd.DataFrame:
    """Copied verbatim from fig2_trait_directness_stats.py."""
    m = members.copy()
    cons = m["consequence"].astype(str)
    m["is_utr"] = cons.str.contains(r"(?:^|,)[35]_prime_UTR_variant(?:$|,)", regex=True)
    m["is_synonymous"] = cons.str.contains(r"(?:^|,)synonymous_variant(?:$|,)", regex=True)
    m["is_noncoding_dna"] = (
        (m["consequence_category"] == "other_noncoding")
        | ((m["consequence_category"] == "synonymous_or_utr")
           & m["is_utr"] & ~m["is_synonymous"])
    )
    nc = m[m["is_noncoding_dna"]].copy()

    def flag(col: str) -> pd.Series:
        return nc[col].astype(str).str.lower() == "true"

    pip = nc["normalized_susie_pip"].to_numpy()
    promoter_resolved = nc["promoter_context_status"] == "resolved"
    lineage_resolved = nc["lineage_accessibility_status"] == "resolved"
    promoter_prox = flag("promoter_proximal")
    abc = flag("abc_enhancer_overlap")
    accessible = flag("lineage_accessible")

    nc["noncoding_pip_mass"] = pip
    nc["promoter_mass"] = pip * (promoter_resolved & promoter_prox)
    nc["abc_mass"] = pip * abc
    nc["accessible_mass"] = pip * (lineage_resolved & accessible)
    nc["resolved_no_context_mass"] = pip * (
        promoter_resolved & lineage_resolved & ~promoter_prox & ~abc & ~accessible)
    nc["unresolved_context_mass"] = pip * (~promoter_resolved | ~lineage_resolved)

    cols = ["noncoding_pip_mass"] + CTX_COLS
    return nc.groupby(["credible_set_uid", "trait_scope", "locus"],
                      as_index=False)[cols].sum()


def cluster_table(frame: pd.DataFrame, num_col: str, den_col: str) -> pd.DataFrame:
    """One row per (scope, operational cluster): mean numerator and denominator.

    Averaging credible-set instances within a cluster is what makes a cluster
    recovered by many studies count once, matching the published estimator.
    """
    g = frame.groupby(["trait_scope", "cluster_id"], as_index=False)[[num_col, den_col]].mean()
    return g.rename(columns={num_col: "num", den_col: "den"})


def analyse(tab: pd.DataFrame, rng: np.random.Generator) -> dict:
    """Fractions, difference, ratio, cluster bootstrap, and a label-shuffle permutation."""
    a = tab[tab["trait_scope"] == SCOPES[0]]
    b = tab[tab["trait_scope"] == SCOPES[1]]
    na, nb = len(a), len(b)
    if na < 3 or nb < 3:
        return {"n_clusters_direct": na, "n_clusters_enzyme": nb}

    def frac(f: pd.DataFrame) -> float:
        d = f["den"].sum()
        return float(f["num"].sum() / d) if d > 0 else np.nan

    fa, fb = frac(a), frac(b)

    an, ad = a["num"].to_numpy(), a["den"].to_numpy()
    bn, bd = b["num"].to_numpy(), b["den"].to_numpy()
    ia = rng.integers(0, na, size=(N_BOOT, na))
    ib = rng.integers(0, nb, size=(N_BOOT, nb))
    da = an[ia].sum(axis=1) / ad[ia].sum(axis=1)
    db = bn[ib].sum(axis=1) / bd[ib].sum(axis=1)
    diff = da - db
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(db > 0, da / db, np.nan)

    # Permutation: shuffle the trait-class label across all scope-cluster rows.
    pool_n = np.concatenate([an, bn])
    pool_d = np.concatenate([ad, bd])
    obs = abs(fa - fb)
    hits = 0
    for _ in range(N_PERM):
        idx = rng.permutation(len(pool_n))
        sa, sb = idx[:na], idx[na:]
        dena, denb = pool_d[sa].sum(), pool_d[sb].sum()
        if dena <= 0 or denb <= 0:
            continue
        hits += abs(pool_n[sa].sum() / dena - pool_n[sb].sum() / denb) >= obs

    return {
        "n_clusters_direct": na, "n_clusters_enzyme": nb,
        "fraction_direct": fa, "fraction_enzyme": fb,
        "difference": fa - fb,
        "diff_boot_lo": float(np.percentile(diff, 2.5)),
        "diff_boot_hi": float(np.percentile(diff, 97.5)),
        "ratio": fa / fb if fb > 0 else np.nan,
        "ratio_boot_lo": float(np.nanpercentile(ratio, 2.5)),
        "ratio_boot_hi": float(np.nanpercentile(ratio, 97.5)),
        # published-style bootstrap p, reproduced for auditability
        "p_bootstrap_diff": float(2 * min((diff <= 0).mean(), (diff >= 0).mean())),
        "p_permutation": (1 + hits) / (1 + N_PERM),
    }


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
    members_path = base / "credible_set_members_annotated.tsv"

    arch = pd.read_csv(arch_path, sep="\t")
    members = pd.read_csv(members_path, sep="\t", low_memory=False)
    arch = arch[arch["trait_scope"].isin(SCOPES)].copy()
    members = members[members["trait_scope"].isin(SCOPES)].copy()

    for (scope, col), expected in PUBLISHED.items():
        got = arch.loc[arch["trait_scope"] == scope, col].mean()
        if not np.isclose(got, expected, rtol=1e-6, atol=1e-9):
            raise AssertionError(f"published {col} for {scope}: expected {expected}, got {got}")
    print("[mult] published instance-weighted architecture reproduced exactly", flush=True)

    arch = assign_clusters(parse_locus(arch))
    arch["coding_mass"] = arch[CODING_COLS].sum(axis=1)
    arch["arch_total"] = arch[ARCH_COLS].sum(axis=1)

    per_set = assign_clusters(parse_locus(build_context_masses(members)))

    specs = [("coding_mass", arch, "coding_mass", "arch_total")] + [
        (c, per_set, c, "noncoding_pip_mass") for c in CTX_COLS]

    rows = []
    for name, frame, num, den in specs:
        res = analyse(cluster_table(frame, num, den), rng)
        res["test"] = name
        res["label"] = TEST_LABEL[name]
        rows.append(res)
        print(f"[mult] {name:<26} direct {res.get('fraction_direct', float('nan')):.4f} "
              f"enzyme {res.get('fraction_enzyme', float('nan')):.4f} "
              f"ratio {res.get('ratio', float('nan')):.3f} "
              f"p_perm {res.get('p_permutation', float('nan')):.4f} "
              f"p_boot {res.get('p_bootstrap_diff', float('nan')):.4g}", flush=True)

    res = pd.DataFrame(rows)
    res["family"] = "fig2_trait_directness_six_tests"
    res["family_size"] = len(res)
    res["q_permutation_BH"] = bh(res["p_permutation"].to_numpy())
    res["q_bootstrap_BH"] = bh(res["p_bootstrap_diff"].to_numpy())
    res["survives_BH_0.05_permutation"] = res["q_permutation_BH"] < 0.05
    res["survives_BH_0.05_bootstrap"] = res["q_bootstrap_BH"] < 0.05

    cols = ["test", "label", "n_clusters_direct", "n_clusters_enzyme",
            "fraction_direct", "fraction_enzyme", "difference", "ratio",
            "ratio_boot_lo", "ratio_boot_hi", "p_permutation", "q_permutation_BH",
            "survives_BH_0.05_permutation", "p_bootstrap_diff", "q_bootstrap_BH",
            "survives_BH_0.05_bootstrap", "family", "family_size"]
    res = res[cols]
    res.to_csv(out / "fig2_family_multiplicity.tsv", sep="\t", index=False)
    print("\n" + res.drop(columns=["family", "family_size", "label"]).to_string(index=False), flush=True)

    manifest = [{"role": r, "path": str(p), "sha256": sha256(p)}
                for r, p in [("architecture", arch_path), ("members", members_path)]]
    (out / "input_manifest.json").write_text(json.dumps(manifest, indent=2))
    (out / "environment.txt").write_text(
        f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
        f"seed {SEED}\nn_boot {N_BOOT}\nn_perm {N_PERM}\n"
        + subprocess.run([sys.executable, "-m", "pip", "freeze"],
                         capture_output=True, text=True).stdout)
    print(f"[mult] wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
