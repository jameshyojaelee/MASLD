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

TWO NULLS, BOTH REPORTED (superseded 2026-09-23; see CORRECTION below).
`fig2_trait_directness_stats.py` reports
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

CORRECTION 2026-09-23 (review items 3 and 4,
`docs/technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23.md`).

Credible-set gate. Only sets with `pip_sum_gate_passed` enter, keyed by
`credible_set_uid` on both the architecture and member tables. The 52 failing
sets (1 direct, 51 enzyme) carry missing mass that pandas summed as 0. A set with
missing or non-finite mass is excluded and written to
`excluded_credible_sets.tsv`, never counted as zero mass.

Physical cluster. Clusters are formed once on the pooled positions of both trait
classes, so a region reached by both classes is one physical cluster holding one
direct row and one enzyme row. After the gate the coding test has 283 rows in
264 physical clusters; 19 are shared, and they hold 19 of the 24 direct clusters.

Bootstrap. Physical clusters are resampled with replacement; a drawn shared
cluster brings both its rows. The percentile interval comes from this bootstrap.
Figure 2 reports intervals only. `nominal_p_bootstrap_diff` and its six-test BH
`nominal_q_bootstrap_BH` are percentile-inversion diagnostics floored at 2/(n+1),
not tests, and no pass/fail column is written.

No permutation p. The candidate null exchanges the class label only among
class-exclusive clusters and keeps the shared pairs fixed. It would relabel 5 of
24 direct clusters while the 19 fixed shared clusters hold 79% of the direct
denominator and 76% of the direct coding mass, so it tests the five exclusive
direct clusters, not the direct-versus-enzyme contrast. Its exchangeability
assumption also fails by ascertainment: enzyme GWAS are larger and reach
smaller-effect loci, so an enzyme-exclusive cluster is not a draw from the same
population as a direct-exclusive one. `p_permutation` is therefore NA and no
permutation q-value is reported.

Read-only inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import zlib
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

# Gate failures in run_19783321, by credible_set_uid (review item 3).
EXPECTED_GATE_FAILURES = {"tier1_direct_masld_pdff": 1, "tier2_liver_enzyme": 51}

# Coding-mass estimate on gated sets, reproduced 2026-09-23 by a standalone
# grouping that does not import this module (ungated it gives 0.22542908161519876
# / 0.034420056715428145 / 6.549352416207798, the superseded published values).
CORRECTED_CODING = {
    "fraction_direct": 0.22680460858995022,
    "fraction_enzyme": 0.03389892569783767,
    "ratio": 6.690613461075479,
}

# (direct rows, enzyme rows, physical clusters, shared clusters) after the gate.
# Coding census reproduced independently; context census from build_context_masses.
EXPECTED_CENSUS = {"coding_mass": (24, 259, 264, 19)}
EXPECTED_CENSUS.update({c: (22, 256, 261, 17) for c in CTX_COLS})

P_PERMUTATION_STATUS = ("not_computed: most direct clusters are shared with the enzyme "
                        "class; no defensible label exchange (see "
                        "fig2_multiplicity_correction.py docstring)")


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


def fit_rng(fit_id: str) -> np.random.Generator:
    """Generator seeded by SEED and a stable fit id, so a fit's draws do not
    depend on which fits ran before it."""
    return np.random.default_rng([SEED, zlib.crc32(fit_id.encode())])


def gate_passed(flag: pd.Series) -> pd.Series:
    """`pip_sum_gate_passed` as bool. pandas reads the TSV's `true`/`false` as
    bool; text spellings such as 'True' or 'TRUE' are accepted too."""
    if flag.dtype == bool:
        return flag
    return flag.astype(str).str.strip().str.lower() == "true"


def apply_gate(arch: pd.DataFrame, members: pd.DataFrame
               ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Keep credible sets that pass the PIP-sum gate and have finite mass.

    The unit is `credible_set_uid`; a set dropped from one table is dropped from
    the other. Returns the two filtered tables and one row per excluded set with
    its reason: `pip_sum_gate_failed`, `nonfinite_architecture_mass`,
    `nonpositive_architecture_mass`, or `nonfinite_member_pip`.
    """
    uid = "credible_set_uid"
    member_flag = gate_passed(members["pip_sum_gate_passed"]).groupby(members[uid]).all()
    arch_flag = gate_passed(arch["pip_sum_gate_passed"]).set_axis(arch[uid])
    if not arch_flag.index.is_unique:
        raise ValueError("credible_set_uid is not unique in the architecture table")
    common = member_flag.index.intersection(arch_flag.index)
    if (member_flag[common] != arch_flag[common]).any():
        raise ValueError("pip_sum_gate_passed disagrees between the two tables")

    failed = set(arch_flag.index[~arch_flag]) | set(member_flag.index[~member_flag])
    passing = arch[~arch[uid].isin(failed)]
    mass = passing[ARCH_COLS].to_numpy(dtype=float)
    nonfinite_arch = set(passing.loc[~np.isfinite(mass).all(axis=1), uid])
    nonpositive_arch = set(passing.loc[np.isfinite(mass).all(axis=1)
                                       & (mass.sum(axis=1) <= 0), uid])
    pip = pd.to_numeric(members["normalized_susie_pip"], errors="coerce")
    bad_pip = ~np.isfinite(pip.to_numpy(dtype=float)) & ~members[uid].isin(failed).to_numpy()
    nonfinite_members = set(members.loc[bad_pip, uid])

    reasons = {}
    for reason, ids in [("nonfinite_member_pip", nonfinite_members),
                        ("nonpositive_architecture_mass", nonpositive_arch),
                        ("nonfinite_architecture_mass", nonfinite_arch),
                        ("pip_sum_gate_failed", failed)]:
        reasons.update(dict.fromkeys(ids, reason))  # later entries win: gate first
    meta = pd.concat([arch[[uid, "trait_scope", "study"]],
                      members[[uid, "trait_scope", "study"]]]).drop_duplicates(uid)
    excluded = meta[meta[uid].isin(reasons)].copy()
    excluded["reason"] = excluded[uid].map(reasons)
    excluded = excluded.sort_values(["reason", "trait_scope", uid]).reset_index(drop=True)
    return (arch[~arch[uid].isin(reasons)].copy(),
            members[~members[uid].isin(reasons)].copy(), excluded)


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


def physical_wide(tab: pd.DataFrame) -> pd.DataFrame:
    """One row per physical cluster with the direct and enzyme num/den side by
    side. A class absent from a cluster contributes 0 to both sums, which leaves
    that class's ratio of sums unchanged."""
    wide = tab.pivot(index="cluster_id", columns="trait_scope", values=["num", "den"])
    wide = wide.reindex(columns=pd.MultiIndex.from_product([["num", "den"], SCOPES]))
    present = wide["den"].notna()
    out = pd.DataFrame({
        "num_a": wide[("num", SCOPES[0])], "den_a": wide[("den", SCOPES[0])],
        "num_b": wide[("num", SCOPES[1])], "den_b": wide[("den", SCOPES[1])],
    }).fillna(0.0)
    out["has_a"] = present[SCOPES[0]].to_numpy()
    out["has_b"] = present[SCOPES[1]].to_numpy()
    return out


def cluster_bootstrap(wide: pd.DataFrame, rng: np.random.Generator,
                      n_boot: int = N_BOOT) -> tuple[np.ndarray, np.ndarray]:
    """Direct and enzyme fractions over `n_boot` resamples of physical clusters.

    One index draw selects whole physical clusters, so a shared cluster's direct
    and enzyme rows always enter or leave a replicate together.
    """
    idx = rng.integers(0, len(wide), size=(n_boot, len(wide)))
    with np.errstate(divide="ignore", invalid="ignore"):
        da = wide["num_a"].to_numpy()[idx].sum(axis=1) / wide["den_a"].to_numpy()[idx].sum(axis=1)
        db = wide["num_b"].to_numpy()[idx].sum(axis=1) / wide["den_b"].to_numpy()[idx].sum(axis=1)
    return da, db


def analyse(tab: pd.DataFrame, rng: np.random.Generator, n_boot: int = N_BOOT) -> dict:
    """Fractions, difference, ratio, and a physical-cluster bootstrap.

    `tab` rows are (trait_scope, cluster_id) with cluster ids formed on the
    pooled positions of both classes, so a shared cluster id is one region. No
    permutation p is computed; see the module docstring.
    """
    a = tab[tab["trait_scope"] == SCOPES[0]]
    b = tab[tab["trait_scope"] == SCOPES[1]]
    na, nb = len(a), len(b)
    if na < 3 or nb < 3:
        return {"n_clusters_direct": na, "n_clusters_enzyme": nb}

    def frac(f: pd.DataFrame) -> float:
        d = f["den"].sum()
        return float(f["num"].sum() / d) if d > 0 else np.nan

    fa, fb = frac(a), frac(b)

    wide = physical_wide(tab)
    da, db = cluster_bootstrap(wide, rng, n_boot)
    ok = np.isfinite(da) & np.isfinite(db)
    diff = (da - db)[ok]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(db[ok] > 0, da[ok] / db[ok], np.nan)
    tail = min(int((diff <= 0).sum()), int((diff >= 0).sum()))

    return {
        "n_clusters_direct": na, "n_clusters_enzyme": nb,
        "n_physical_clusters": len(wide),
        "n_shared_clusters": int((wide["has_a"] & wide["has_b"]).sum()),
        "fraction_direct": fa, "fraction_enzyme": fb,
        "difference": fa - fb,
        "diff_boot_lo": float(np.percentile(diff, 2.5)),
        "diff_boot_hi": float(np.percentile(diff, 97.5)),
        "ratio": fa / fb if fb > 0 else np.nan,
        "ratio_boot_lo": float(np.nanpercentile(ratio, 2.5)),
        "ratio_boot_hi": float(np.nanpercentile(ratio, 97.5)),
        "n_boot_finite": int(ok.sum()),
        # Percentile-inversion diagnostic from the physical-cluster bootstrap, not
        # a test: twice the smaller tail share of the difference at zero, floored
        # at 2/(n+1) so a tail with no draws is not reported as 0. Figure 2
        # reports intervals only and makes no significance call.
        "nominal_p_bootstrap_diff": min(1.0, 2 * (1 + tail) / (1 + len(diff))),
        "p_permutation": np.nan,
        "p_permutation_status": P_PERMUTATION_STATUS,
    }


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

    arch, members, excluded = apply_gate(arch, members)
    excluded.to_csv(out / "excluded_credible_sets.tsv", sep="\t", index=False)
    gate_fail = (excluded[excluded["reason"] == "pip_sum_gate_failed"]
                 .groupby("trait_scope").size().to_dict())
    if gate_fail != EXPECTED_GATE_FAILURES:
        raise AssertionError(f"gate failures {gate_fail}, expected {EXPECTED_GATE_FAILURES}")
    print(f"[mult] excluded credible sets by reason: "
          f"{excluded.groupby(['reason', 'trait_scope']).size().to_dict()}", flush=True)

    arch = assign_clusters(parse_locus(arch))
    arch["coding_mass"] = arch[CODING_COLS].sum(axis=1)
    arch["arch_total"] = arch[ARCH_COLS].sum(axis=1)

    per_set = assign_clusters(parse_locus(build_context_masses(members)))

    specs = [("coding_mass", arch, "coding_mass", "arch_total")] + [
        (c, per_set, c, "noncoding_pip_mass") for c in CTX_COLS]

    rows = []
    for name, frame, num, den in specs:
        res = analyse(cluster_table(frame, num, den), fit_rng(name))
        res["test"] = name
        res["label"] = TEST_LABEL[name]
        rows.append(res)
        census = (res["n_clusters_direct"], res["n_clusters_enzyme"],
                  res["n_physical_clusters"], res["n_shared_clusters"])
        if census != EXPECTED_CENSUS[name]:
            raise AssertionError(f"{name} census (direct, enzyme, physical, shared) "
                                 f"{census}, expected {EXPECTED_CENSUS[name]}")
        print(f"[mult] {name:<26} direct {res.get('fraction_direct', float('nan')):.4f} "
              f"enzyme {res.get('fraction_enzyme', float('nan')):.4f} "
              f"ratio {res.get('ratio', float('nan')):.3f} "
              f"[{res['ratio_boot_lo']:.3f}, {res['ratio_boot_hi']:.3f}] "
              f"physical {census[2]} shared {census[3]} "
              f"nominal_p_boot {res.get('nominal_p_bootstrap_diff', float('nan')):.4g}", flush=True)

    coding = rows[0]
    for key, expected in CORRECTED_CODING.items():
        tol = 1e-9 if key.startswith("fraction") else 1e-7 * expected
        if abs(coding[key] - expected) > tol:
            raise AssertionError(f"coding {key}: expected {expected}, got {coding[key]}")
    print("[mult] gated coding estimate matches the independent reproduction", flush=True)

    res = pd.DataFrame(rows)
    res["family"] = "fig2_trait_directness_six_tests"
    res["family_size"] = len(res)
    res["inference_unit"] = "physical_1Mb_cluster_pooled_over_trait_classes"
    res["q_permutation_BH"] = np.nan
    # BH over the six tests applied to the nominal diagnostic; no pass/fail flag
    res["nominal_q_bootstrap_BH"] = bh(res["nominal_p_bootstrap_diff"].to_numpy())
    res["seed"] = SEED
    res["n_boot"] = N_BOOT

    cols = ["test", "label", "n_clusters_direct", "n_clusters_enzyme",
            "n_physical_clusters", "n_shared_clusters",
            "fraction_direct", "fraction_enzyme", "difference",
            "diff_boot_lo", "diff_boot_hi", "ratio",
            "ratio_boot_lo", "ratio_boot_hi", "n_boot_finite",
            "nominal_p_bootstrap_diff", "nominal_q_bootstrap_BH",
            "p_permutation", "q_permutation_BH", "p_permutation_status",
            "inference_unit", "family", "family_size", "seed", "n_boot"]
    res = res[cols]
    res.to_csv(out / "fig2_family_multiplicity.tsv", sep="\t", index=False)
    print("\n" + res.drop(columns=["family", "family_size", "label", "p_permutation_status",
                                   "inference_unit"]).to_string(index=False), flush=True)

    manifest = [{"role": r, "path": str(p), "sha256": sha256(p)}
                for r, p in [("architecture", arch_path), ("members", members_path)]]
    (out / "input_manifest.json").write_text(json.dumps(manifest, indent=2))
    (out / "environment.txt").write_text(
        f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
        f"seed {SEED} (per test: fit_rng(test))\nn_boot {N_BOOT}\nn_perm not used\n"
        + subprocess.run([sys.executable, "-m", "pip", "freeze"],
                         capture_output=True, text=True).stdout)
    print(f"[mult] wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
