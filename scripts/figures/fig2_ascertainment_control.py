#!/usr/bin/env python
"""Ascertainment control for the direct-trait coding-concentration claim.

THE CLAIM UNDER TEST. Fine-mapped PIP mass at direct hepatic-fat loci is
coding-concentrated relative to liver-enzyme loci: 22.4% vs 3.4% of annotated
mass, a ratio of 6.6x on independent loci (`fig2_trait_directness_stats.py`).

THE OBJECTION. Direct-trait GWAS are less powered than liver-enzyme GWAS. An
underpowered study only reaches significance for large-effect variants, and
large-effect variants are disproportionately protein-altering because coding
variation escapes the regulatory buffering that constrains common noncoding
effects. On that account the 6.6x measures study power, not liver biology.

The objection is not idle here, and the direction is worse than it looks from
raw sample sizes. The direct-trait arm contains the portfolio's largest studies
by N_tot (2021_34841290_NAFLD_EUR at 778,614), but most are case-control designs
with few cases, so their EFFECTIVE sample size is small: 2023_36280732_NAFLD_UKBB
has 651 cases in 397,020 participants, an effective N near 2,600. The liver
enzyme arm is quantitative, where effective N is the full sample: UKBB_ALT is
343,850. That is a two-order-of-magnitude power gap pointing exactly the way the
objection requires.

THE CONTROL. For each credible set, take its lead variant's standardized effect
size from the study's own summary statistics,

    beta_std  =  |z| / sqrt(N_eff),        z = beta / se

which is scale-free across binary and quantitative traits and separates the size
of an effect from the power to detect it. Then ask whether the trait-class
difference survives when direct and enzyme credible sets are compared at matched
beta_std, and whether the trait-class term survives adjustment for it.

N_eff is N_tot for quantitative traits and the standard 4 / (1/n_cases +
1/n_controls) for binary ones.

Read-only inputs. Every quantity is estimated with the independent locus as the
unit, matching the claim it is testing, with a cluster bootstrap over loci.
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
LOCUS_GAP_BP = 1_000_000
CALIPER_SD = 0.2
CHUNK = 2_000_000

SCOPES = ["tier1_direct_masld_pdff", "tier2_liver_enzyme"]
SCOPE_LABEL = {
    "tier1_direct_masld_pdff": "Direct MASLD / liver fat",
    "tier2_liver_enzyme": "Liver enzymes",
}
CODING_COLS = ["protein_altering_pip_mass", "canonical_splice_pip_mass"]
ARCH_COLS = CODING_COLS + ["synonymous_or_utr_pip_mass", "other_noncoding_pip_mass"]

# Published instance-weighted values from the read-only candidate. Asserted, not
# assumed, so this script cannot silently drift off the claim it is testing.
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


def parse_locus(frame: pd.DataFrame) -> pd.DataFrame:
    """Copied verbatim from fig2_trait_directness_stats.py."""
    text = frame["locus"].astype(str)
    frame = frame.copy()
    frame["chrom"] = text.str.split(".").str[0]
    frame["pos"] = pd.to_numeric(text.str.split(".").str[1], errors="coerce")
    if frame["pos"].isna().any():
        raise ValueError("unparsable locus positions present")
    return frame


def assign_independent_loci(frame: pd.DataFrame) -> pd.DataFrame:
    """Copied verbatim from fig2_trait_directness_stats.py.

    Single-linkage clustering of credible-set positions within a chromosome, run
    across the whole portfolio so one physical locus recovered by several studies
    collapses to one unit.
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


def effective_n(row: pd.Series) -> float:
    """N_eff: N_tot for quantitative traits, 4/(1/cases + 1/controls) for binary.

    The binary form is the standard effective-sample-size correction. It is the
    reason this control is necessary: a 397,020-participant study with 651 cases
    has the power of roughly 2,600 balanced observations, not 397,020.
    """
    n_tot = float(row["N_tot"])
    if str(row["trait_type"]).strip().lower() == "binary":
        n_cases = float(row["N_cases"])
        n_ctrl = n_tot - n_cases
        if n_cases <= 0 or n_ctrl <= 0:
            return np.nan
        return 4.0 / (1.0 / n_cases + 1.0 / n_ctrl)
    return n_tot


def lead_variants(members: pd.DataFrame) -> pd.DataFrame:
    """One row per credible set: the maximum-PIP member."""
    m = members.sort_values("normalized_susie_pip", ascending=False)
    lead = m.drop_duplicates("credible_set_uid", keep="first")
    return lead[["credible_set_uid", "study", "trait_scope", "locus",
                 "chromosome", "position", "allele1", "allele2",
                 "normalized_susie_pip"]].copy()


def fetch_effects(lead: pd.DataFrame, registry: pd.DataFrame,
                  sumstats_dir: Path) -> pd.DataFrame:
    """Stream each study's summary statistics and pull its lead variants' beta/se.

    Matching is on (chromosome, position) with an allele check in both
    orientations, because the fine-mapping export and the reformatted summary
    statistics do not guarantee a common allele order. A lead variant that cannot
    be matched is returned with NaN and counted, never silently dropped.
    """
    out = []
    for study, block in lead.groupby("study", sort=True):
        reg = registry[registry["study_name"] == study]
        if reg.empty:
            print(f"[effects] {study}: absent from registry, skipped", flush=True)
            continue
        path = sumstats_dir / Path(str(reg.iloc[0]["sumstats_path"])).name
        if not path.exists():
            print(f"[effects] {study}: {path.name} not on disk, skipped", flush=True)
            continue

        want = block.copy()
        want["chromosome"] = want["chromosome"].astype(str).str.replace("^chr", "", regex=True)
        # String keys with a vectorised `isin` rather than a Python membership
        # loop: these files total 29 GB, so a per-row loop over every chunk would
        # dominate the runtime.
        keys = set(want["chromosome"] + ":" + want["position"].astype("int64").astype(str))
        found = {}
        for chunk in pd.read_csv(
                path, sep="\t", chunksize=CHUNK, low_memory=False,
                usecols=["chromosome", "position", "allele1", "allele2", "beta", "se"]):
            chunk["chromosome"] = chunk["chromosome"].astype(str).str.replace("^chr", "", regex=True)
            chunk = chunk[pd.to_numeric(chunk["position"], errors="coerce").notna()]
            chunk["position"] = chunk["position"].astype("int64")
            ck = chunk["chromosome"] + ":" + chunk["position"].astype(str)
            hit = chunk[ck.isin(keys)]
            for r in hit.itertuples(index=False):
                found.setdefault((r.chromosome, r.position), []).append(
                    (str(r.allele1), str(r.allele2), r.beta, r.se))

        betas, ses, matched = [], [], []
        for r in want.itertuples(index=False):
            cands = found.get((r.chromosome, int(r.position)), [])
            pick = None
            for a1, a2, b, s in cands:
                if {a1.upper(), a2.upper()} == {str(r.allele1).upper(), str(r.allele2).upper()}:
                    pick = (b, s)
                    break
            if pick is None and len(cands) == 1:
                # position-unique match with an allele mismatch: keep it but flag
                pick = (cands[0][2], cands[0][3])
                matched.append("position_only")
            elif pick is None:
                matched.append("unmatched")
            else:
                matched.append("exact")
            betas.append(np.nan if pick is None else pick[0])
            ses.append(np.nan if pick is None else pick[1])
        want["beta"] = betas
        want["se"] = ses
        want["match_type"] = matched
        n_ok = int(np.isfinite(want["beta"]).sum())
        print(f"[effects] {study}: {n_ok}/{len(want)} lead variants matched", flush=True)
        out.append(want)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def locus_means(frame: pd.DataFrame, value: str) -> pd.DataFrame:
    """Credible sets -> loci (mean), so a locus counts once regardless of how
    many studies recovered it. This is the unit the 6.6x claim uses."""
    return frame.groupby(["trait_scope", "locus_id"], as_index=False)[value].mean()


def ratio_with_ci(frame: pd.DataFrame, value: str, rng: np.random.Generator) -> dict:
    """Direct/enzyme ratio of locus-mean coding fraction, cluster bootstrap over loci."""
    lm = locus_means(frame, value)
    a = lm.loc[lm["trait_scope"] == SCOPES[0], value].to_numpy()
    b = lm.loc[lm["trait_scope"] == SCOPES[1], value].to_numpy()
    if len(a) < 3 or len(b) < 3:
        return {"n_loci_direct": len(a), "n_loci_enzyme": len(b)}
    point = a.mean() / b.mean() if b.mean() > 0 else np.nan
    draws = np.empty(N_BOOT)
    for i in range(N_BOOT):
        aa = a[rng.integers(0, len(a), len(a))]
        bb = b[rng.integers(0, len(b), len(b))]
        draws[i] = aa.mean() / bb.mean() if bb.mean() > 0 else np.nan
    draws = draws[np.isfinite(draws)]
    # Two-sided permutation p for the ratio, labels shuffled over loci.
    pool = np.concatenate([a, b])
    obs = abs(np.log(point)) if np.isfinite(point) and point > 0 else np.nan
    perm = 0
    for _ in range(N_BOOT):
        idx = rng.permutation(len(pool))
        pa, pb = pool[idx[:len(a)]], pool[idx[len(a):]]
        if pb.mean() > 0 and pa.mean() > 0:
            perm += abs(np.log(pa.mean() / pb.mean())) >= obs
    return {
        "n_loci_direct": len(a), "n_loci_enzyme": len(b),
        "mean_direct": a.mean(), "mean_enzyme": b.mean(),
        "ratio": point,
        "ci_lo": float(np.quantile(draws, 0.025)),
        "ci_hi": float(np.quantile(draws, 0.975)),
        "perm_p_two_sided": (1 + perm) / (1 + N_BOOT),
    }


def caliper_match(frame: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """1:1 nearest-neighbour match of enzyme to direct credible sets on log10 beta_std.

    Matching is at credible-set level because beta_std is a property of a
    study-by-locus observation, then the matched set is collapsed to loci for
    estimation, so the unit of inference is still the independent locus.
    """
    d = frame.dropna(subset=["log_beta_std"]).copy()
    cases = d[d["trait_scope"] == SCOPES[0]]
    pool = d[d["trait_scope"] == SCOPES[1]]
    if cases.empty or pool.empty:
        return pd.DataFrame()
    caliper = CALIPER_SD * d["log_beta_std"].std(ddof=1)
    pool_sorted = pool.sort_values("log_beta_std")
    pv = pool_sorted["log_beta_std"].to_numpy()
    pidx = pool_sorted.index.to_numpy()
    used: set = set()
    rows = []
    for i, row in cases.sample(frac=1.0, random_state=int(rng.integers(1 << 31))).iterrows():
        order = np.argsort(np.abs(pv - row["log_beta_std"]))
        for k in order:
            if pidx[k] in used:
                continue
            if abs(pv[k] - row["log_beta_std"]) > caliper:
                break
            used.add(pidx[k])
            rows.append(i)
            rows.append(pidx[k])
            break
    return d.loc[rows]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    root = Path(args.project_root)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    base = root / ("GWAS/finemapping/results/candidates/"
                   "noncoding-dna-precoloc-2026-08-12/run_19783321")
    arch_path = base / "credible_set_pip_architecture.tsv"
    members_path = base / "credible_set_members_annotated.tsv"
    registry_path = root / "GWAS/finemapping/config/gwas_registry.tsv"
    sumstats_dir = root / "GWAS/finemapping/data/sumstats"

    arch = pd.read_csv(arch_path, sep="\t")
    members = pd.read_csv(members_path, sep="\t", low_memory=False)
    registry = pd.read_csv(registry_path, sep="\t")
    arch = arch[arch["trait_scope"].isin(SCOPES)].copy()
    members = members[members["trait_scope"].isin(SCOPES)].copy()

    # ---- reproduce the published instance-weighted values as an assertion ----
    for (scope, col), expected in PUBLISHED.items():
        got = arch.loc[arch["trait_scope"] == scope, col].mean()
        if not np.isclose(got, expected, rtol=1e-6, atol=1e-9):
            raise AssertionError(
                f"published instance-weighted {col} for {scope}: expected {expected}, got {got}")
    print("[control] published instance-weighted architecture reproduced exactly", flush=True)

    arch["coding_mass"] = arch[CODING_COLS].sum(axis=1)
    arch = assign_independent_loci(parse_locus(arch))

    # ---- baseline: the claim being tested, on independent loci ---------------
    baseline = ratio_with_ci(arch, "coding_mass", rng)
    print(f"[control] BASELINE locus-weighted coding ratio = {baseline['ratio']:.3f} "
          f"({baseline['n_loci_direct']} direct loci, {baseline['n_loci_enzyme']} enzyme)",
          flush=True)

    # ---- effect sizes and effective N ---------------------------------------
    registry["n_eff"] = registry.apply(effective_n, axis=1)
    lead = lead_variants(members)
    eff = fetch_effects(lead, registry, sumstats_dir)
    eff = eff.merge(registry[["study_name", "n_eff", "N_tot", "N_cases", "trait_type"]],
                    left_on="study", right_on="study_name", how="left")
    eff["z"] = eff["beta"] / eff["se"]
    eff["abs_z"] = eff["z"].abs()
    eff["beta_std"] = eff["abs_z"] / np.sqrt(eff["n_eff"])
    eff["log_beta_std"] = np.log10(eff["beta_std"].where(eff["beta_std"] > 0))
    eff.to_csv(out / "lead_variant_effects.tsv", sep="\t", index=False)

    merged = arch.merge(
        eff[["credible_set_uid", "beta", "se", "z", "abs_z", "n_eff",
             "beta_std", "log_beta_std", "match_type"]],
        on="credible_set_uid", how="left")
    merged.to_csv(out / "credible_sets_with_power.tsv", sep="\t", index=False)

    cov = merged.groupby("trait_scope")["log_beta_std"].apply(lambda s: s.notna().mean())
    print("[control] lead-variant effect coverage by scope:\n" + cov.to_string(), flush=True)

    # ---- the imbalance the objection depends on ------------------------------
    imb = merged.groupby("trait_scope").agg(
        n_credible_sets=("credible_set_uid", "size"),
        median_n_eff=("n_eff", "median"),
        median_abs_z=("abs_z", "median"),
        median_beta_std=("beta_std", "median"),
        q25_beta_std=("beta_std", lambda s: s.quantile(0.25)),
        q75_beta_std=("beta_std", lambda s: s.quantile(0.75)),
    ).reset_index()
    imb["scope_label"] = imb["trait_scope"].map(SCOPE_LABEL)
    imb.to_csv(out / "power_imbalance.tsv", sep="\t", index=False)
    print("[control] power imbalance:\n" + imb.to_string(index=False), flush=True)

    # ---- control 1: matched on standardized effect size ---------------------
    matched = caliper_match(merged, rng)
    matched_res = ratio_with_ci(matched, "coding_mass", rng) if not matched.empty else {}
    if not matched.empty:
        smd = ((matched.loc[matched.trait_scope == SCOPES[0], "log_beta_std"].mean()
                - matched.loc[matched.trait_scope == SCOPES[1], "log_beta_std"].mean())
               / matched["log_beta_std"].std(ddof=1))
        matched_res["match_smd_log_beta_std"] = smd
        matched_res["n_matched_credible_sets"] = len(matched)
        print(f"[control] MATCHED ratio = {matched_res.get('ratio', float('nan')):.3f} "
              f"(SMD after matching {smd:+.3f})", flush=True)

    # ---- control 2: adjustment, cluster-robust over loci --------------------
    reg_rows = []
    try:
        import statsmodels.api as sm
        d = merged.dropna(subset=["log_beta_std", "coding_mass"]).copy()
        d["is_direct"] = (d["trait_scope"] == SCOPES[0]).astype(float)
        d["log_n_eff"] = np.log10(d["n_eff"])
        for name, cols in (("unadjusted", ["is_direct"]),
                           ("adj_beta_std", ["is_direct", "log_beta_std"]),
                           ("adj_beta_std_and_n", ["is_direct", "log_beta_std", "log_n_eff"])):
            X = sm.add_constant(d[cols], has_constant="add")
            fit = sm.OLS(d["coding_mass"], X).fit(
                cov_type="cluster", cov_kwds={"groups": d["locus_id"]})
            reg_rows.append({
                "model": name,
                "beta_is_direct": fit.params["is_direct"],
                "se_is_direct": fit.bse["is_direct"],
                "p_is_direct": fit.pvalues["is_direct"],
                "ci_lo": fit.conf_int().loc["is_direct", 0],
                "ci_hi": fit.conf_int().loc["is_direct", 1],
                "n_obs": int(fit.nobs), "n_clusters": d["locus_id"].nunique(),
            })
    except ImportError:
        print("[control] statsmodels unavailable; regression arm skipped", flush=True)
    if reg_rows:
        pd.DataFrame(reg_rows).to_csv(out / "adjustment_models.tsv", sep="\t", index=False)
        print("[control] adjustment models:\n"
              + pd.DataFrame(reg_rows).to_string(index=False), flush=True)

    summary = pd.DataFrame([
        {"analysis": "baseline_locus_weighted", **baseline},
        {"analysis": "matched_on_standardized_effect", **matched_res},
    ])
    summary.to_csv(out / "ascertainment_control_summary.tsv", sep="\t", index=False)

    manifest = [{"role": r, "path": str(p), "sha256": sha256(p)}
                for r, p in [("architecture", arch_path), ("members", members_path),
                             ("gwas_registry", registry_path)]]
    (out / "input_manifest.json").write_text(json.dumps(manifest, indent=2))
    (out / "environment.txt").write_text(
        f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
        f"seed {SEED}\nn_boot {N_BOOT}\ncaliper_sd {CALIPER_SD}\n"
        + subprocess.run([sys.executable, "-m", "pip", "freeze"],
                         capture_output=True, text=True).stdout)
    print(f"[control] wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
