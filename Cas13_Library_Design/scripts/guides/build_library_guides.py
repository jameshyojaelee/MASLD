#!/usr/bin/env python
"""
build_library_guides.py — assemble the full Cas13 sgRNA guide library.

Reads the v8 target roster (cas13_library.csv), pulls N guides/gene for the
protein-coding + lncRNA targets from the cached parquet index (no miRNA tier in
v8), merges the library annotations, and pulls control-gene guides (essential
genes; positive controls are in-library targets, not re-pulled).

Outputs (Cas13_Library_Design/data/guides/):
  cas13_library_guides_<release>.csv   target guides + annotations
  cas13_control_guides_<release>.csv   essential-gene QC guides
  cas13_guides_coverage_report.csv           per-gene coverage (PC + lncRNA roster)
  GUIDE_BUILD_MANIFEST.txt                   provenance + parameters + counts

Non-targeting + safe-harbor controls are NOT generated here (non-genic; not in
the upstream pool) — see GUIDE_DESIGN.md "Controls TODO".

Run on a compute node (queries are light; the index build is the heavy step).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guides_config as cfg
import guide_selection as gs


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(cfg.PROJECT), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# Essentiality (DepMap Chronos, liver cell lines) — annotated onto the guide table
# ---------------------------------------------------------------------------
def liver_essentiality() -> pd.DataFrame:
    """Per-human-gene mean DepMap Chronos across LIVER cell lines (Model.csv
    OncotreeLineage == 'Liver'); mirrors RNA-seq/27a_assemble_evidence_atlas.R layer 7.

    Returns: gene_symbol_human, essentiality_chronos_liver, n_liver_lines,
    is_essential_liver (Chronos < cfg.ESSENTIAL_CHRONOS_THR; more negative = more
    essential). The 412 MB gene-effect matrix is read only on a cache miss; the small
    per-gene table is cached to cfg.LIVER_ESSENTIALITY_CACHE thereafter.
    """
    cache = cfg.LIVER_ESSENTIALITY_CACHE
    if cache.exists():
        return pd.read_csv(cache)
    models = pd.read_csv(cfg.DEPMAP_MODEL, usecols=["ModelID", "OncotreeLineage"])
    liver_ids = set(models.loc[models["OncotreeLineage"] == "Liver", "ModelID"])
    eff = pd.read_csv(cfg.DEPMAP_GENE_EFFECT, index_col=0)      # rows=ModelID, cols 'SYM (ENTREZ)'
    eff = eff.loc[eff.index.isin(liver_ids)]
    chronos = eff.mean(axis=0, skipna=True)                    # per-gene mean across liver lines
    n_lines = eff.notna().sum(axis=0).astype(int)
    syms = chronos.index.to_series().str.replace(r"\s*\(\d+\)$", "", regex=True)
    tbl = (pd.DataFrame({
                "gene_symbol_human": syms.values,
                "essentiality_chronos_liver": chronos.round(4).values,
                "n_liver_lines": n_lines.values,
            })
           .dropna(subset=["gene_symbol_human"])
           .query("gene_symbol_human != ''")
           .drop_duplicates("gene_symbol_human"))
    tbl["is_essential_liver"] = tbl["essentiality_chronos_liver"] < cfg.ESSENTIAL_CHRONOS_THR
    cache.parent.mkdir(parents=True, exist_ok=True)
    tbl.to_csv(cache, index=False)
    print(f"[build] liver essentiality: {len(eff)} liver lines x {len(tbl)} genes "
          f"-> {cache.name} (cached)")
    return tbl


# ---------------------------------------------------------------------------
# Targets (protein-coding + lncRNA; miRNA deferred)
# ---------------------------------------------------------------------------
def build_targets(lazy, lib: pd.DataFrame, n: int):
    tgt = lib[lib["biotype"].isin(cfg.TARGET_BIOTYPES)].copy()
    targets = pd.DataFrame({
        "query": tgt["gene_symbol_mouse"].values,
        "gene_id_base": tgt["gene_id_mouse"].values,
    })
    guides, summary = gs.select_guides(lazy, targets, n=n)

    annot = lib.set_index("gene_id_mouse")[cfg.LIBRARY_ANNOT_COLS]
    guides = guides.merge(annot, left_on="gene_id_mouse", right_index=True, how="left")
    summary = summary.merge(annot, left_on="gene_id_mouse", right_index=True, how="left")
    return guides, summary


# ---------------------------------------------------------------------------
# Controls (positive controls + essential genes)
# ---------------------------------------------------------------------------
def build_controls(lazy, n: int):
    frames = []
    notes = {}

    # v7: positive controls are now folded into the TARGET roster (flagged
    # is_positive_control in cas13_library.csv), so they are NOT pulled
    # separately here -- that would double-list them. Only the assay-QC controls
    # (essential genes) are generated; non-targeting + safe-harbor remain a TODO.
    notes["positive_controls"] = "in-library targets (is_positive_control flag); not re-pulled"

    # Essential genes (Cas13-activity QC)
    es_map = gs.human_to_mouse_ids(cfg.ESSENTIAL_GENES_HUMAN)
    es_guides, es_summary = gs.select_guides(
        lazy, es_map[["query", "gene_id_base"]], n=n)
    if not es_guides.empty:
        id2human = dict(zip(es_map["gene_id_base"], es_map["query"]))
        es_guides["human_symbol"] = es_guides["gene_id_mouse"].map(id2human)
        es_guides["expected_direction"] = "depletion"
        es_guides["control_type"] = "essential"
        frames.append(es_guides)
    notes["essential_input"] = len(cfg.ESSENTIAL_GENES_HUMAN)
    notes["essential_mapped"] = int(es_map["mapped"].sum())
    notes["essential_with_guides"] = int((es_summary["n_selected"] > 0).sum())

    controls = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return controls, notes


# ---------------------------------------------------------------------------
# Coverage report: complete accounting of all library genes (PC + lncRNA; no miRNA in v8)
# ---------------------------------------------------------------------------
def _load_posthoc(release: str) -> pd.DataFrame:
    """Per-gene provenance from the upstream funnel (keyed by unversioned id).

    `constitutive_library_status` + `genome_n_guides_passed` are the trustworthy
    'why' columns; `exclusion_reason` is a noisy per-transcript-stage annotation
    kept only for reference (see GUIDE_DESIGN.md).
    """
    ph_path = cfg.RELEASES[release]["posthoc"]
    want = ["gene_id", "constitutive_library_status",
            "genome_n_guides_passed", "exclusion_reason"]
    ph = pd.read_csv(ph_path, sep="\t", usecols=want, dtype=str)
    ph["gene_id_base"] = ph["gene_id"].str.split(".").str[0]
    return ph.drop(columns=["gene_id"]).drop_duplicates("gene_id_base")


def coverage_report(lib: pd.DataFrame, summary: pd.DataFrame, n: int,
                    release: str) -> pd.DataFrame:
    # Start from the library roster so identity (symbol/biotype) is always present,
    # even for genes absent from the upstream pool.
    target = lib[lib["biotype"].isin(cfg.TARGET_BIOTYPES)][
        ["gene_id_mouse", "gene_symbol_mouse", "biotype"]].copy()
    stats = summary[["gene_id_mouse", "n_available_pool", "n_selected",
                     "short_flag", "missing_flag"]]
    cov = target.merge(stats, on="gene_id_mouse", how="left")
    cov["n_available_pool"] = cov["n_available_pool"].fillna(0).astype(int)
    cov["n_selected"] = cov["n_selected"].fillna(0).astype(int)
    cov["missing_flag"] = cov["missing_flag"].fillna(True).astype(bool)
    cov["short_flag"] = cov["short_flag"].fillna(True).astype(bool)
    cov["status"] = "built"
    cov.loc[cov["missing_flag"], "status"] = "no_guides_in_source"
    cov.loc[cov["short_flag"] & ~cov["missing_flag"], "status"] = f"short_lt_{n}"

    # Merge upstream funnel provenance (why a gene is guide-poor / absent).
    cov = cov.merge(_load_posthoc(release), left_on="gene_id_mouse",
                    right_on="gene_id_base", how="left").drop(columns=["gene_id_base"])

    # miRNA targets are deferred (not in the upstream PCG+lncRNA pool).
    mir = lib[lib["biotype"] == "miRNA"][["gene_id_mouse", "gene_symbol_mouse", "biotype"]].copy()
    mir["n_available_pool"] = 0
    mir["n_selected"] = 0
    mir["short_flag"] = True
    mir["missing_flag"] = True
    mir["status"] = "deferred_miRNA"
    return pd.concat([cov, mir], ignore_index=True)


def write_manifest(path: Path, release: str, n: int,
                   lib: pd.DataFrame, target_guides: pd.DataFrame,
                   target_summary: pd.DataFrame, control_notes: dict,
                   n_control_guides: int):
    bt = lib["biotype"].value_counts().to_dict()
    by_tier = (target_guides.drop_duplicates("gene_id_mouse")["tier"].value_counts().to_dict()
               if "tier" in target_guides.columns else {})
    n_short = int((target_summary["short_flag"] & ~target_summary["missing_flag"]).sum())
    n_missing = int(target_summary["missing_flag"].sum())
    lines = [
        "Cas13 sgRNA Guide Library — BUILD MANIFEST",
        "=" * 50,
        f"built_utc           : {datetime.now(timezone.utc).isoformat()}",
        f"git_sha             : {_git_sha()}",
        f"release             : {release}",
        f"source_genome_csv   : {cfg.RELEASES[release]['genome']}",
        f"source_posthoc_tsv  : {cfg.RELEASES[release]['posthoc']}",
        f"index_parquet       : {cfg.index_path(release)}",
        f"library_roster      : {cfg.LIBRARY_CSV} ({cfg.LIBRARY_VERSION})",
        "",
        "PARAMETERS",
        f"  guides_per_gene   : {n}",
        f"  guide_len         : {cfg.GUIDE_LEN}",
        f"  score_guard       : TIGER>={cfg.TIGER_MIN} OR Cas13Design>={cfg.CAS13_MIN}",
        f"  homopolymer_guard : {', '.join(cfg.HOMOPOLYMERS)}",
        f"  ortholog_min_tier : {cfg.ORTHOLOG_MIN_TIER}",
        f"  selection         : constitutive-first (max isoform coverage), transcript-spread, combined_score",
        "",
        "LIBRARY ROSTER (cas13_library.csv)",
        f"  total_genes       : {len(lib)}",
        *[f"    {k:16s}: {v}" for k, v in bt.items()],
        "",
        "TARGET GUIDES (protein_coding + lncRNA; miRNA deferred)",
        f"  target_genes      : {target_summary['biotype'].isin(cfg.TARGET_BIOTYPES).sum()}",
        f"  guides_written    : {len(target_guides)}",
        f"  genes_short_<{n}    : {n_short}",
        f"  genes_no_guides   : {n_missing}",
        *[f"    tier {k:14s}: {v} genes" for k, v in by_tier.items()],
        "",
        "CONTROL GUIDES",
        *[f"  {k:28s}: {v}" for k, v in control_notes.items()],
        f"  control_guides_written      : {n_control_guides}",
        "",
        "DEFERRED / TODO",
        f"  miRNA_targets_deferred      : {int((lib['biotype']=='miRNA').sum())} (need separate pri-miRNA design)",
        "  non_targeting_controls      : TODO (scrambled; not in upstream pool)",
        "  safe_harbor_controls        : TODO (Rosa26/AAVS1-equiv; not in upstream pool)",
        "",
    ]
    path.write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", default=cfg.DEFAULT_RELEASE,
                    choices=list(cfg.RELEASES.keys()))
    ap.add_argument("--n", type=int, default=cfg.N_GUIDES_DEFAULT)
    args = ap.parse_args()

    cfg.GUIDES_DIR.mkdir(parents=True, exist_ok=True)
    lazy = gs.load_index(args.release)
    lib = pd.read_csv(cfg.LIBRARY_CSV)
    print(f"[build] library roster: {len(lib)} genes "
          f"({lib['biotype'].value_counts().to_dict()})")

    # --- liver essentiality (DepMap Chronos) annotation, joined on gene_symbol_human ---
    ess = liver_essentiality()

    # --- targets ---
    target_guides, target_summary = build_targets(lazy, lib, args.n)
    target_guides = target_guides.merge(ess, on="gene_symbol_human", how="left")
    n_ess_genes = target_guides.dropna(subset=["essentiality_chronos_liver"])["gene_id_mouse"].nunique()
    out_targets = cfg.GUIDES_DIR / f"cas13_library_guides_{args.release}.csv"
    target_guides.to_csv(out_targets, index=False)
    print(f"[build] targets: {len(target_guides)} guides over "
          f"{target_guides['gene_id_mouse'].nunique()} genes -> {out_targets.name}"
          f"  (liver-essentiality scored: {n_ess_genes} genes)")

    # --- controls ---
    controls, control_notes = build_controls(lazy, args.n)
    if not controls.empty and "human_symbol" in controls.columns:
        controls = controls.merge(
            ess.rename(columns={"gene_symbol_human": "human_symbol"}),
            on="human_symbol", how="left")
    out_controls = cfg.GUIDES_DIR / f"cas13_control_guides_{args.release}.csv"
    controls.to_csv(out_controls, index=False)
    print(f"[build] controls: {len(controls)} guides -> {out_controls.name}  {control_notes}")

    # --- coverage report (full PC + lncRNA roster; no miRNA in v8) ---
    cov = coverage_report(lib, target_summary, args.n, args.release)
    out_cov = cfg.GUIDES_DIR / "cas13_guides_coverage_report.csv"
    cov.to_csv(out_cov, index=False)
    print(f"[build] coverage: {cov['status'].value_counts().to_dict()} -> {out_cov.name}")

    # --- persistent unguideable list (genes with NO qualifying guide in the pool) ---
    # rebuild_cas13_library.R reads this to drop untargetable genes at the SOURCE.
    # Accumulated as a UNION across runs so a gene already dropped from the roster is
    # not "forgotten" when a later run no longer processes it (convergence guarantee).
    ung_now = cov.loc[cov["status"] == "no_guides_in_source",
                      ["gene_id_mouse", "gene_symbol_mouse", "biotype",
                       "genome_n_guides_passed", "exclusion_reason"]].copy()
    ung_path = cfg.GUIDES_DIR / f"unguideable_{args.release}.csv"
    if ung_path.exists():
        prev = pd.read_csv(ung_path, dtype=str)
        ung = pd.concat([prev, ung_now.astype(str)], ignore_index=True) \
                .drop_duplicates("gene_id_mouse", keep="last")
    else:
        ung = ung_now
    ung.to_csv(ung_path, index=False)
    print(f"[build] unguideable: {len(ung_now)} this run, {len(ung)} cumulative -> {ung_path.name}")

    # --- manifest ---
    write_manifest(cfg.GUIDES_DIR / "GUIDE_BUILD_MANIFEST.txt", args.release, args.n,
                   lib, target_guides, target_summary, control_notes, len(controls))

    # --- invariants ---
    bad_len = target_guides.loc[target_guides["guide_seq"].str.len() != cfg.GUIDE_LEN]
    dup_id = target_guides["guide_id"].duplicated().sum()
    print(f"[build] INVARIANTS: non-{cfg.GUIDE_LEN}bp guides={len(bad_len)} | "
          f"duplicate guide_ids={dup_id}")
    if len(bad_len):
        print(f"[build] WARNING: {len(bad_len)} guides not {cfg.GUIDE_LEN}bp", file=sys.stderr)
    print("[build] DONE")


if __name__ == "__main__":
    main()
