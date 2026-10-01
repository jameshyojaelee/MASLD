#!/usr/bin/env python3
"""Recompute the genetics x disease-state interface on the PROMOTED COLOC release.

Why this exists. Two headline claims rested on the retired July colocalization
set (473 SuSiE genes, 447 jointly testable, 14,931-gene background):

  independence   threshold-free Spearman rho = 0.019 (p = 0.232) between
                 colocalization strength and disease-state effect
  no buffering   expression-matched McNemar ladder is flat; 23 colocalized genes
                 strongly DE against 24 matched controls at the 0.50 floor

Both are functions of the colocalized gene set AND of the jointly testable
background, and the background is DEFINED as the intersection of the bulk-tested
and colocalization-tested universes. A new COLOC release therefore moves the
denominator, not just the numerator. That is exactly the asymmetric-denominator
trap that made an earlier version of this analysis wrong, so nothing is carried
over from July: the background is rebuilt from the promoted release.

Two bulk substrates are run side by side and never pooled:

  official    27,638-gene read-count universe (the substrate July used)
  fragment    23,370-gene five-cohort fragment universe (validated, and per
              RESULTS.md the intended substrate, but NOT promoted at time of
              writing)

Reporting both is deliberate. If the conclusion holds on both, it does not
depend on an unpromoted substrate choice; if it differs, that difference is
itself the result and must be seen before either number enters prose.
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
from scipy import stats

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
COLOC = PROJ / "GWAS/finemapping/results/susie_coloc/gene_level_coloc_tier12.csv"
BULK = {
    "official_27638": PROJ / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
                           / "results/integration/canonical_deg_results.csv",
    "fragment_23370": PROJ / "RNA-seq/results/manuscript_release/candidates"
                           / "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION"
                           / "arms/F_five/results/integration/deg_results.csv",
}

SEED = 20260817
CALIPER_SD = 0.2
N_BOOT = 2000
FLOORS = [0.00, 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.00]

# Canonical DEG gate, per STATUS.md's 2026-08-12 migration.
CANON_PADJ = 0.05
CANON_LFC = 0.50


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def base_id(s: pd.Series) -> pd.Series:
    return s.astype(str).str.split(".").str[0]


def load_coloc(path: Path = COLOC) -> pd.DataFrame:
    d = pd.read_csv(path)
    d["g"] = base_id(d["ensembl"])
    # One duplicated ensembl exists; keep the strongest evidence per gene so the
    # gene set is a set, not a bag.
    d = (d.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False)
           .drop_duplicates("g", keep="first"))
    return d


def build(bulk_path: Path, coloc: pd.DataFrame) -> dict:
    deg = pd.read_csv(bulk_path)
    deg["g"] = base_id(deg["gene"])
    deg = deg.drop_duplicates("g", keep="first")

    # Joint testability rebuilt from the promoted release: a gene is jointly
    # testable when it was assayed in the bulk contrast AND appears in the
    # colocalization universe. No July flag is reused.
    # Genes with NA best SuSiE PP.H4 stay in this universe. A corrected release marks
    # a tier-1/2 SuSiE test that failed coloc's shared-posterior check as
    # susie_state_t12 == "untestable"; such a gene counts as not SuSiE-positive, the
    # same as a tested gene at PP.H4 <= 0.5, because only former positives were
    # rechecked.
    keep = ["g", "coloc_best_susie_pp4", "coloc_best_abf_pp4"]
    keep += [c for c in ("susie_state_t12",) if c in coloc.columns]
    jt = deg.merge(coloc[keep], on="g", how="inner")
    jt["susie"] = pd.to_numeric(jt["coloc_best_susie_pp4"], errors="coerce")
    jt["abf"] = pd.to_numeric(jt["coloc_best_abf_pp4"], errors="coerce")
    jt["is_deg"] = (jt["padj"] < CANON_PADJ) & (jt["logFC"].abs() > CANON_LFC)
    jt["abs_lfc"] = jt["logFC"].abs()
    jt["abs_shrunk"] = jt["shrunk_logFC"].abs()
    return {"deg": deg, "jt": jt}


def overlap_stats(jt: pd.DataFrame, genetic: pd.Series, label: str) -> dict:
    n = len(jt)
    a = int((genetic & jt["is_deg"]).sum())
    b = int((genetic & ~jt["is_deg"]).sum())
    c = int((~genetic & jt["is_deg"]).sum())
    d = int((~genetic & ~jt["is_deg"]).sum())
    exp = (a + b) * (a + c) / n if n else np.nan
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return {
        "genetic_definition": label,
        "n_joint_testable": n,
        "n_genetic_joint": a + b,
        "n_deg_joint": a + c,
        "n_overlap_observed": a,
        "n_overlap_expected": exp,
        "pct_genetic_not_deg": 100 * b / (a + b) if (a + b) else np.nan,
        "fisher_or": orr,
        "fisher_p": p,
    }


def threshold_free(jt: pd.DataFrame, col: str) -> dict:
    m = jt[col].notna() & jt["shrunk_logFC"].notna()
    x = jt.loc[m, col].to_numpy()
    signed = jt.loc[m, "shrunk_logFC"].to_numpy()
    absol = np.abs(signed)
    rs, ps = stats.spearmanr(x, signed)
    ra, pa = stats.spearmanr(x, absol)
    m_raw = jt[col].notna() & jt["logFC"].notna()
    x_raw = jt.loc[m_raw, col].to_numpy()
    signed_raw = jt.loc[m_raw, "logFC"].to_numpy()
    absol_raw = np.abs(signed_raw)
    rs_raw, ps_raw = stats.spearmanr(x_raw, signed_raw)
    ra_raw, pa_raw = stats.spearmanr(x_raw, absol_raw)
    return {"pp4_column": col, "n_pairs": int(m.sum()),
            "spearman_signed": rs, "p_signed": ps,
            "spearman_absolute": ra, "p_absolute": pa,
            "n_pairs_raw": int(m_raw.sum()),
            "spearman_signed_raw": rs_raw, "p_signed_raw": ps_raw,
            "spearman_absolute_raw": ra_raw, "p_absolute_raw": pa_raw}


def matched_ladder(
    jt: pd.DataFrame, genetic: pd.Series, rng: np.random.Generator
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """Expression-matched McNemar ladder.

    The unmatched ladder declines with the |log2FC| floor because colocalized
    genes are higher-expressed and a high floor preferentially admits
    low-expression genes whose effects are inflated conditional on passing FDR.
    Matching on AveExpr removes that, and the ladder is then scored WITHIN pairs
    so the discordant-pair test is the estimand the ladder actually makes.
    """
    d = jt.dropna(subset=["AveExpr", "logFC", "padj", "treat_fdr"]).copy()
    cases = d[genetic.reindex(d.index).fillna(False)]
    pool = d[~genetic.reindex(d.index).fillna(False)]
    if len(cases) < 20 or len(pool) < len(cases):
        return pd.DataFrame(), {}, pd.DataFrame()

    caliper = CALIPER_SD * d["AveExpr"].std(ddof=1)
    used: set[int] = set()
    pairs = []
    pool_sorted = pool.sort_values("AveExpr")
    pv = pool_sorted["AveExpr"].to_numpy()
    pidx = pool_sorted.index.to_numpy()
    for i, row in cases.sample(frac=1.0, random_state=int(rng.integers(1 << 31))).iterrows():
        j = np.searchsorted(pv, row["AveExpr"])
        order = sorted(range(len(pv)), key=lambda k: abs(pv[k] - row["AveExpr"]))
        for k in order:
            cand = pidx[k]
            if cand in used:
                continue
            if abs(pv[k] - row["AveExpr"]) > caliper:
                break
            used.add(cand)
            pairs.append((i, cand))
            break
    if not pairs:
        return pd.DataFrame(), {}, pd.DataFrame()

    ci = [p[0] for p in pairs]
    co = [p[1] for p in pairs]
    case_d = d.loc[ci]
    ctrl_d = d.loc[co]
    smd = ((case_d["AveExpr"].mean() - ctrl_d["AveExpr"].mean())
           / np.sqrt((case_d["AveExpr"].var(ddof=1) + ctrl_d["AveExpr"].var(ddof=1)) / 2))

    rows = []
    for f in FLOORS:
        ca = ((case_d["padj"] < CANON_PADJ) & (case_d["logFC"].abs() > f)).to_numpy()
        co_ = ((ctrl_d["padj"] < CANON_PADJ) & (ctrl_d["logFC"].abs() > f)).to_numpy()
        b = int((ca & ~co_).sum())
        c = int((~ca & co_).sum())
        orr = np.nan if c == 0 else b / c
        p = stats.binomtest(b, b + c, 0.5).pvalue if (b + c) else np.nan
        rows.append({"floor": f, "n_pairs": len(pairs), "match_smd_aveexpr": smd,
                     "n_case_de": int(ca.sum()), "n_control_de": int(co_.sum()),
                     "discordant_case_only": b, "discordant_control_only": c,
                     "mcnemar_or": orr, "mcnemar_p": p})
    case_treat = (case_d["treat_fdr"] < 0.05).to_numpy()
    ctrl_treat = (ctrl_d["treat_fdr"] < 0.05).to_numpy()
    case_only = int((case_treat & ~ctrl_treat).sum())
    control_only = int((~case_treat & ctrl_treat).sum())
    both = int((case_treat & ctrl_treat).sum())
    neither = int((~case_treat & ~ctrl_treat).sum())
    treat_or = np.nan if control_only == 0 else case_only / control_only
    treat_p = (stats.binomtest(case_only, case_only + control_only, 0.5).pvalue
               if (case_only + control_only) else np.nan)
    treat_summary = {
        "n_pairs": len(pairs),
        "match_smd_aveexpr": smd,
        "treat_lfc": float(case_d["treat_lfc"].iloc[0]),
        "n_case_treat": int(case_treat.sum()),
        "n_control_treat": int(ctrl_treat.sum()),
        "discordant_case_only": case_only,
        "discordant_control_only": control_only,
        "both_treat": both,
        "neither_treat": neither,
        "mcnemar_or": treat_or,
        "mcnemar_p": treat_p,
    }
    pair_table = pd.DataFrame({
        "genetic_gene": case_d["g"].to_numpy(),
        "matched_control_gene": ctrl_d["g"].to_numpy(),
        "genetic_AveExpr": case_d["AveExpr"].to_numpy(),
        "matched_control_AveExpr": ctrl_d["AveExpr"].to_numpy(),
        "genetic_logFC": case_d["logFC"].to_numpy(),
        "matched_control_logFC": ctrl_d["logFC"].to_numpy(),
        "genetic_treat_fdr": case_d["treat_fdr"].to_numpy(),
        "matched_control_treat_fdr": ctrl_d["treat_fdr"].to_numpy(),
        "genetic_treat_supported": case_treat,
        "matched_control_treat_supported": ctrl_treat,
    })
    return pd.DataFrame(rows), treat_summary, pair_table


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--coloc", type=Path, default=COLOC,
                    help="gene_level_coloc_tier12.csv of the COLOC release (default: adopted)")
    ap.add_argument("--bulk-fragment", type=Path, default=BULK["fragment_23370"],
                    help="deg_results.csv of the five-cohort fragment fit (default: F_five v6)")
    args = ap.parse_args()
    # The official arm stays in so the fragment arm's matching uses the same
    # draws of the seeded generator as the adopted run.
    bulk_arms = dict(BULK, fragment_23370=args.bulk_fragment)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    coloc = load_coloc(args.coloc)
    manifest = [{"role": "coloc_tier12", "path": str(args.coloc), "sha256": sha256(args.coloc),
                 "n_genes": int(coloc["g"].nunique())}]

    ov_rows, tf_rows, ladders, treat_rows, pair_tables = [], [], [], [], []
    for arm, path in bulk_arms.items():
        if not path.exists():
            print(f"SKIP {arm}: {path} absent", file=sys.stderr)
            continue
        manifest.append({"role": f"bulk_{arm}", "path": str(path), "sha256": sha256(path)})
        b = build(path, coloc)
        jt = b["jt"]
        print(f"\n=== {arm}: bulk universe {len(b['deg'])}, jointly testable {len(jt)} ===")

        defs = {
            "susie_pp4_0.5": jt["susie"] >= 0.5,
            "susie_or_abf_pp4_0.5": (jt["susie"] >= 0.5) | (jt["abf"] >= 0.5),
        }
        for label, mask in defs.items():
            mask = mask.fillna(False)
            r = overlap_stats(jt, mask, label)
            r["bulk_arm"] = arm
            r["n_genes_bulk"] = len(b["deg"])   # the arm label is an identifier, not a count
            if "susie_state_t12" in jt:
                # Reported as their own state; they are in n_joint_testable and not genetic.
                r["n_susie_untestable_joint"] = int((jt["susie_state_t12"] == "untestable").sum())
            ov_rows.append(r)
            print(f"  {label}: genetic {r['n_genetic_joint']}, DEG {r['n_deg_joint']}, "
                  f"overlap {r['n_overlap_observed']} vs expected {r['n_overlap_expected']:.1f}, "
                  f"OR {r['fisher_or']:.3f} p {r['fisher_p']:.3g}")

            lad, treat_summary, pair_table = matched_ladder(jt, mask, rng)
            if not lad.empty:
                lad["bulk_arm"] = arm
                lad["genetic_definition"] = label
                ladders.append(lad)
                treat_summary["bulk_arm"] = arm
                treat_summary["n_genes_bulk"] = len(b["deg"])
                treat_summary["genetic_definition"] = label
                treat_rows.append(treat_summary)
                pair_table["bulk_arm"] = arm
                pair_table["genetic_definition"] = label
                pair_tables.append(pair_table)

        for col in ("susie", "abf"):
            t = threshold_free(jt, col)
            t["bulk_arm"] = arm
            t["n_genes_bulk"] = len(b["deg"])
            if col == "susie" and "susie_state_t12" in jt:
                # No PP.H4 to rank, so these genes are outside the threshold-free
                # correlation only; they remain in every overlap denominator.
                t["n_untestable_without_pp4"] = int(
                    ((jt["susie_state_t12"] == "untestable") & jt["susie"].isna()).sum())
            tf_rows.append(t)
            print(f"  threshold-free {col}: rho_signed {t['spearman_signed']:.4f} "
                  f"(p {t['p_signed']:.3g}), rho_abs {t['spearman_absolute']:.4f} "
                  f"(p {t['p_absolute']:.3g}), n {t['n_pairs']}; "
                  f"raw rho_signed {t['spearman_signed_raw']:.4f} "
                  f"(p {t['p_signed_raw']:.3g})")

    pd.DataFrame(ov_rows).to_csv(out / "overlap_by_arm.tsv", sep="\t", index=False)
    pd.DataFrame(tf_rows).to_csv(out / "threshold_free_by_arm.tsv", sep="\t", index=False)
    if ladders:
        pd.concat(ladders).to_csv(out / "matched_ladder_by_arm.tsv", sep="\t", index=False)
        pd.DataFrame(treat_rows).to_csv(
            out / "matched_treat_by_arm.tsv", sep="\t", index=False)
        pd.concat(pair_tables).to_csv(
            out / "matched_pairs_by_arm.tsv", sep="\t", index=False)
    (out / "input_manifest.json").write_text(json.dumps(manifest, indent=2))
    (out / "environment.txt").write_text(
        f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
        f"seed {SEED}\n" +
        subprocess.run([sys.executable, "-m", "pip", "freeze"],
                       capture_output=True, text=True).stdout)
    print(f"\nwrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
