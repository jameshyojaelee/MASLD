#!/usr/bin/env python
"""Region-balanced re-estimation of fine-mapped PIP architecture by phenotype directness.

The published candidate reports credible-set-instance-weighted means. Instances are
study-by-region units, so a region discovered by many studies is counted many times.
This script reproduces the published instance-weighted values as an assertion, then
re-estimates every quantity with operational 1-Mb clusters as the unit and a cluster
bootstrap over them for uncertainty.

NAMING (2026-08-17). These groups are OPERATIONAL CLUSTERS, not independent loci. No
LD is consulted in forming them, a 1-Mb window can span more than one independent
signal, and long-range LD can link clusters. The function name
`assign_independent_loci` and the `locus_id` column are retained so existing outputs
stay joinable, but no output of this script may be described as an independent locus.

MULTIPLICITY (2026-08-17). The six direct-vs-enzyme comparisons this script feeds
(coding mass plus five noncoding context categories) are one family and are corrected
together in `fig2_multiplicity_correction.py`. Note also that `p_two_sided_diff_gt0`
below is a BOOTSTRAP p, asking whether a difference clears zero, not whether the
trait-class grouping produces it. It is not the right test for a group contrast on its
own; see `GWAS/finemapping/results/fig2_multiplicity/20260817T161128Z/MULTIPLICITY.md`.

Read-only inputs. Writes source tables for the Figure 2 companion panels.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260813
N_BOOT = 10000
LOCUS_GAP_BP = 1_000_000

SCOPES = ["tier1_direct_masld_pdff", "tier2_liver_enzyme"]
SCOPE_LABEL = {
    "tier1_direct_masld_pdff": "Direct MASLD / liver fat",
    "tier2_liver_enzyme": "Liver enzymes",
}

ARCH_COLS = [
    "protein_altering_pip_mass",
    "canonical_splice_pip_mass",
    "synonymous_or_utr_pip_mass",
    "other_noncoding_pip_mass",
]
# Context is defined on NONCODING variants only and normalised by noncoding PIP
# mass, exactly as `scripts/figures/noncoding_dna_plot_upgrade.R` does. The
# per-credible-set `credible_set_regulatory_context.tsv` uses total mass as its
# denominator and is therefore a different estimand; it is not used here.
CTX_EXCLUSIVE = ["resolved_no_context_mass", "unresolved_context_mass"]
CTX_NONEXCLUSIVE = ["promoter_mass", "abc_mass", "accessible_mass"]
CTX_ALL = CTX_NONEXCLUSIVE + CTX_EXCLUSIVE

PUBLISHED_CONTEXT = {
    ("tier1_direct_masld_pdff", "promoter_mass"): 0.0787249975787247,
    ("tier2_liver_enzyme", "promoter_mass"): 0.0445912158230019,
    ("tier1_direct_masld_pdff", "abc_mass"): 0.00230953599846294,
    ("tier2_liver_enzyme", "abc_mass"): 0.0156986213596511,
    ("tier1_direct_masld_pdff", "accessible_mass"): 0.0731733768551714,
    ("tier2_liver_enzyme", "accessible_mass"): 0.0415012498200219,
    ("tier1_direct_masld_pdff", "resolved_no_context_mass"): 0.334108741121616,
    ("tier2_liver_enzyme", "resolved_no_context_mass"): 0.306842957487844,
    ("tier1_direct_masld_pdff", "unresolved_context_mass"): 0.521129683271921,
    ("tier2_liver_enzyme", "unresolved_context_mass"): 0.604308028252816,
}

# Published instance-weighted values in the immutable candidate; asserted, not assumed.
PUBLISHED = {
    ("tier1_direct_masld_pdff", "protein_altering_pip_mass"): 0.157736519057125,
    ("tier2_liver_enzyme", "protein_altering_pip_mass"): 0.0325683421565398,
    ("tier1_direct_masld_pdff", "canonical_splice_pip_mass"): 0.000361875951271015,
    ("tier2_liver_enzyme", "canonical_splice_pip_mass"): 0.0,
    ("tier1_direct_masld_pdff", "synonymous_or_utr_pip_mass"): 0.0845868982467048,
    ("tier2_liver_enzyme", "synonymous_or_utr_pip_mass"): 0.0511001327057017,
    ("tier1_direct_masld_pdff", "other_noncoding_pip_mass"): 0.757314706744875,
    ("tier2_liver_enzyme", "other_noncoding_pip_mass"): 0.916331525137759,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_locus(frame: pd.DataFrame) -> pd.DataFrame:
    """Split the stored `locus` field into chromosome and base-pair position."""
    text = frame["locus"].astype(str)
    frame = frame.copy()
    frame["chrom"] = text.str.split(".").str[0]
    frame["pos"] = pd.to_numeric(text.str.split(".").str[1], errors="coerce")
    if frame["pos"].isna().any():
        raise ValueError("unparsable locus positions present")
    return frame


def assign_independent_loci(frame: pd.DataFrame) -> pd.DataFrame:
    """Single-linkage clustering of credible-set positions within a chromosome.

    Two credible sets join the same locus when their sentinel positions are within
    LOCUS_GAP_BP. Clustering is done across the whole portfolio so that the same
    physical locus recovered by several studies collapses to one unit.
    """
    frame = frame.sort_values(["chrom", "pos"]).copy()
    keys = []
    for chrom, block in frame.groupby("chrom", sort=False):
        pos = block["pos"].to_numpy()
        new_cluster = np.concatenate([[True], np.diff(pos) > LOCUS_GAP_BP])
        keys.append(pd.Series(
            [f"{chrom}:{i}" for i in np.cumsum(new_cluster)], index=block.index))
    frame["locus_id"] = pd.concat(keys).reindex(frame.index)
    return frame


def build_context_masses(members: pd.DataFrame) -> pd.DataFrame:
    """Per-credible-set noncoding context masses, reproducing the v5 definitions.

    `is_noncoding_dna` and every context predicate are copied verbatim from
    scripts/figures/noncoding_dna_plot_upgrade.R lines 58-172.
    """
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

    cols = ["noncoding_pip_mass"] + CTX_ALL
    per_set = nc.groupby(
        ["credible_set_uid", "trait_scope", "locus"], as_index=False)[cols].sum()
    return per_set


def context_locus_estimates(per_set: pd.DataFrame, cols: list[str],
                            rng: np.random.Generator) -> pd.DataFrame:
    """Locus-weighted context fractions with a cluster bootstrap over loci.

    Each operational 1-Mb region contributes the mean of its credible-set instances,
    so a region recovered by many studies counts once. Fractions are ratios of summed
    region-level masses, matching the published sum/sum estimator within each unit.
    """
    per_locus = per_set.groupby("locus_id")[cols + ["noncoding_pip_mass"]].mean()
    num = per_locus[cols].to_numpy()
    den = per_locus["noncoding_pip_mass"].to_numpy()
    n = num.shape[0]
    idx = rng.integers(0, n, size=(N_BOOT, n))
    draws_num = num[idx].sum(axis=1)
    draws_den = den[idx].sum(axis=1)[:, None]
    fracs = draws_num / draws_den
    out = []
    for j, col in enumerate(cols):
        lo, hi = np.percentile(fracs[:, j], [2.5, 97.5])
        out.append({
            "category": col,
            "instance_weighted_fraction": float(
                per_set[col].sum() / per_set["noncoding_pip_mass"].sum()),
            "locus_weighted_fraction": float(
                per_locus[col].sum() / per_locus["noncoding_pip_mass"].sum()),
            "locus_boot_lo": float(lo),
            "locus_boot_hi": float(hi),
            "absolute_pip_mass_units": float(per_set[col].sum()),
            "n_credible_set_instances": int(len(per_set)),
            "n_operational_regions_1mb": int(n),
            "n_loci_with_any_mass": int((per_locus[col] > 0).sum()),
        })
    return pd.DataFrame(out)


def context_ratio_bootstrap(a: pd.DataFrame, b: pd.DataFrame, col: str,
                            rng: np.random.Generator) -> dict[str, float]:
    """Between-trait ratio and difference for a context fraction, over loci."""
    def draws(frame: pd.DataFrame) -> np.ndarray:
        per_locus = frame.groupby("locus_id")[[col, "noncoding_pip_mass"]].mean()
        num = per_locus[col].to_numpy()
        den = per_locus["noncoding_pip_mass"].to_numpy()
        idx = rng.integers(0, num.shape[0], size=(N_BOOT, num.shape[0]))
        return num[idx].sum(axis=1) / den[idx].sum(axis=1)

    fa, fb = draws(a), draws(b)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(fb > 0, fa / fb, np.nan)
    diff = fa - fb
    return {
        "ratio_median": float(np.nanmedian(ratio)),
        "ratio_lo": float(np.nanpercentile(ratio, 2.5)),
        "ratio_hi": float(np.nanpercentile(ratio, 97.5)),
        "diff_median": float(np.median(diff)),
        "diff_lo": float(np.percentile(diff, 2.5)),
        "diff_hi": float(np.percentile(diff, 97.5)),
        "p_two_sided_diff_gt0": float(2 * min((diff <= 0).mean(), (diff >= 0).mean())),
    }


def instance_weighted(frame: pd.DataFrame, cols: list[str]) -> dict[str, float]:
    """Mass-weighted pooled fraction, matching the published estimator."""
    total = frame[cols].to_numpy().sum()
    return {c: float(frame[c].sum() / total) for c in cols}


def locus_weighted_point(frame: pd.DataFrame, cols: list[str]) -> dict[str, float]:
    """Each operational 1-Mb region contributes once, after averaging its instances."""
    per_locus = frame.groupby("locus_id")[cols].mean()
    total = per_locus.to_numpy().sum()
    return {c: float(per_locus[c].sum() / total) for c in cols}


def locus_bootstrap(frame: pd.DataFrame, cols: list[str], rng: np.random.Generator
                    ) -> dict[str, tuple[float, float]]:
    """Cluster bootstrap resampling operational clusters with replacement."""
    per_locus = frame.groupby("locus_id")[cols].mean()
    mat = per_locus.to_numpy()
    n = mat.shape[0]
    idx = rng.integers(0, n, size=(N_BOOT, n))
    draws = mat[idx].sum(axis=1)
    fracs = draws / draws.sum(axis=1, keepdims=True)
    out = {}
    for j, c in enumerate(cols):
        lo, hi = np.percentile(fracs[:, j], [2.5, 97.5])
        out[c] = (float(lo), float(hi))
    return out


def ratio_bootstrap(a: pd.DataFrame, b: pd.DataFrame, col: str, cols: list[str],
                    rng: np.random.Generator) -> dict[str, float]:
    """Bootstrap the between-trait ratio and difference on the locus-weighted scale."""
    pa = a.groupby("locus_id")[cols].mean().to_numpy()
    pb = b.groupby("locus_id")[cols].mean().to_numpy()
    j = cols.index(col)
    ia = rng.integers(0, pa.shape[0], size=(N_BOOT, pa.shape[0]))
    ib = rng.integers(0, pb.shape[0], size=(N_BOOT, pb.shape[0]))
    da = pa[ia].sum(axis=1)
    db = pb[ib].sum(axis=1)
    fa = da[:, j] / da.sum(axis=1)
    fb = db[:, j] / db.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(fb > 0, fa / fb, np.nan)
    diff = fa - fb
    return {
        "ratio_median": float(np.nanmedian(ratio)),
        "ratio_lo": float(np.nanpercentile(ratio, 2.5)),
        "ratio_hi": float(np.nanpercentile(ratio, 97.5)),
        "diff_median": float(np.median(diff)),
        "diff_lo": float(np.percentile(diff, 2.5)),
        "diff_hi": float(np.percentile(diff, 97.5)),
        "p_two_sided_diff_gt0": float(2 * min((diff <= 0).mean(), (diff >= 0).mean())),
    }


def concentration(frame: pd.DataFrame, col: str, cols: list[str]) -> pd.DataFrame:
    """How few loci carry the mass in a category, ranked."""
    per_locus = frame.groupby("locus_id")[cols].mean()
    share = per_locus[col] / per_locus[col].sum()
    ranked = share.sort_values(ascending=False)
    pos = frame.groupby("locus_id")["pos"].median()
    chrom = frame.groupby("locus_id")["chrom"].first()
    out = pd.DataFrame({
        "locus_id": ranked.index,
        "chromosome": chrom.reindex(ranked.index).to_numpy(),
        "median_position_hg19": pos.reindex(ranked.index).to_numpy(),
        "share_of_category_mass": ranked.to_numpy(),
        "cumulative_share": np.cumsum(ranked.to_numpy()),
    })
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    base = root / ("GWAS/finemapping/results/candidates/"
                   "noncoding-dna-precoloc-2026-08-12/run_19783321")
    arch_path = base / "credible_set_pip_architecture.tsv"
    members_path = base / "credible_set_members_annotated.tsv"

    arch = pd.read_csv(arch_path, sep="\t")
    members = pd.read_csv(members_path, sep="\t", low_memory=False)

    manifest = pd.DataFrame([
        {"input_role": "credible_set_pip_architecture", "path": str(arch_path),
         "sha256": sha256(arch_path), "rows": len(arch)},
        {"input_role": "credible_set_members_annotated", "path": str(members_path),
         "sha256": sha256(members_path), "rows": len(members)},
    ])
    manifest.to_csv(out / "input_manifest.tsv", sep="\t", index=False)

    # Reliability gate exactly as the published candidate applied it.
    arch = arch[arch["pip_sum_gate_passed"].astype(str).str.lower() == "true"]
    members = members[members["pip_sum_gate_passed"].astype(str).str.lower() == "true"]

    ctx = build_context_masses(members)

    # Cluster loci once on the architecture table, then map the same locus ids onto
    # the context table so both families use an identical independence structure.
    arch = assign_independent_loci(parse_locus(arch))
    ctx = ctx.merge(arch[["credible_set_uid", "locus_id"]], on="credible_set_uid",
                    how="left", validate="one_to_one")
    if ctx["locus_id"].isna().any():
        raise ValueError("context rows without a matching architecture locus")

    rng = np.random.default_rng(SEED)
    rows = []
    design_rows = []

    for scope in SCOPES:
        a = arch[arch["trait_scope"] == scope]
        c = ctx[ctx["trait_scope"] == scope]
        design_rows.append({
            "trait_scope": scope,
            "trait_scope_label": SCOPE_LABEL[scope],
            "n_credible_set_instances": len(a),
            "n_operational_regions_1mb": a["locus_id"].nunique(),
            "n_studies": a["study"].nunique(),
            "n_traits": a["trait"].nunique(),
            "instances_per_locus": len(a) / a["locus_id"].nunique(),
            "largest_single_study_share": (
                a["study"].value_counts().iloc[0] / len(a)),
            "largest_single_study": a["study"].value_counts().index[0],
        })

        iw = instance_weighted(a, ARCH_COLS)
        lw = locus_weighted_point(a, ARCH_COLS)
        ci = locus_bootstrap(a, ARCH_COLS, rng)
        for col in ARCH_COLS:
            rows.append({
                "trait_scope": scope,
                "trait_scope_label": SCOPE_LABEL[scope],
                "family": "consequence",
                "denominator": "total_credible_set_pip_mass",
                "category": col,
                "instance_weighted_fraction": iw[col],
                "locus_weighted_fraction": lw[col],
                "locus_boot_lo": ci[col][0],
                "locus_boot_hi": ci[col][1],
                "absolute_pip_mass_units": float(a[col].sum()),
                "n_credible_set_instances": len(a),
                "n_operational_regions_1mb": a["locus_id"].nunique(),
                "n_loci_with_any_mass": int(
                    (a.groupby("locus_id")[col].sum() > 0).sum()),
            })

        ctx_est = context_locus_estimates(c, CTX_ALL, rng)
        ctx_est.insert(0, "trait_scope", scope)
        ctx_est.insert(1, "trait_scope_label", SCOPE_LABEL[scope])
        ctx_est.insert(2, "family", np.where(
            ctx_est["category"].isin(CTX_NONEXCLUSIVE),
            "context_nonexclusive", "context_exclusive"))
        ctx_est.insert(3, "denominator", "noncoding_pip_mass")
        rows.extend(ctx_est.to_dict("records"))

    estimates = pd.DataFrame(rows)
    design = pd.DataFrame(design_rows)

    # Assertion: the instance-weighted consequence values must reproduce the
    # immutable published candidate exactly.
    for (scope, col), expected in PUBLISHED.items():
        got = estimates.loc[
            (estimates.trait_scope == scope) & (estimates.category == col),
            "instance_weighted_fraction"].iloc[0]
        if not np.isclose(got, expected, atol=1e-9, rtol=0):
            raise AssertionError(
                f"instance-weighted {scope}/{col} = {got!r}, published {expected!r}")
    for (scope, col), expected in PUBLISHED_CONTEXT.items():
        got = estimates.loc[
            (estimates.trait_scope == scope) & (estimates.category == col),
            "instance_weighted_fraction"].iloc[0]
        if not np.isclose(got, expected, atol=1e-9, rtol=0):
            raise AssertionError(
                f"instance-weighted context {scope}/{col} = {got!r}, "
                f"published {expected!r}")
    print("PASS: reproduced all published instance-weighted consequence "
          "and regulatory-context values")

    a1 = arch[arch.trait_scope == SCOPES[0]]
    a2 = arch[arch.trait_scope == SCOPES[1]]
    contrasts = [{"family": "consequence", "category": "protein_altering_pip_mass",
                  **ratio_bootstrap(a1, a2, "protein_altering_pip_mass", ARCH_COLS, rng)}]
    c1 = ctx[ctx.trait_scope == SCOPES[0]]
    c2 = ctx[ctx.trait_scope == SCOPES[1]]
    for col in CTX_ALL:
        contrasts.append({"family": "context", "category": col,
                          **context_ratio_bootstrap(c1, c2, col, rng)})
    contrasts = pd.DataFrame(contrasts)

    conc = concentration(a1, "protein_altering_pip_mass", ARCH_COLS)
    conc.insert(0, "trait_scope", SCOPES[0])
    conc.insert(1, "category", "protein_altering_pip_mass")

    estimates.to_csv(out / "pip_architecture_estimates.tsv", sep="\t", index=False)
    design.to_csv(out / "portfolio_independence.tsv", sep="\t", index=False)
    contrasts.to_csv(out / "trait_class_contrasts.tsv", sep="\t", index=False)
    conc.to_csv(out / "direct_trait_coding_concentration.tsv", sep="\t", index=False)

    with open(out / "environment.txt", "w") as handle:
        handle.write(f"python {sys.version}\n")
        handle.write(f"numpy {np.__version__}\npandas {pd.__version__}\n")
        handle.write(f"seed {SEED}\nn_boot {N_BOOT}\nlocus_gap_bp {LOCUS_GAP_BP}\n")
        handle.write(subprocess.run(
            [sys.executable, "-m", "pip", "freeze"], capture_output=True,
            text=True).stdout)

    pd.set_option("display.width", 200)
    print("\n=== portfolio independence ===")
    print(design.to_string(index=False))
    print("\n=== estimates ===")
    print(estimates.round(5).to_string(index=False))
    print("\n=== trait-class contrasts (locus-weighted) ===")
    print(contrasts.round(4).to_string(index=False))
    print("\n=== direct-trait coding mass concentration (top 12 loci) ===")
    print(conc.head(12).round(4).to_string(index=False))
    print(f"\nwrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
