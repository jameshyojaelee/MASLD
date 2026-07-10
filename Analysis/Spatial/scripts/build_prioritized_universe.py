#!/usr/bin/env python
"""Rebuild the prioritized MASLD-target universe (Fig-4 validation substrate).

WHAT THIS IS
------------
Regenerable, auditable builder for the three on-disk universe files under
`Analysis/Spatial/results/universe_validation/`:

    universe_genetic.txt          genetic-evidence arm   (locked count: 3,038)
    universe_transcriptomic.txt   transcriptomic arm     (locked count: 8,260)
    prioritized_universe_FINAL.txt = genetic ∪ transcriptomic   (10,044)

Consumers: `15g_prioritized_gsmap_enrichment.py`, `15h_arm_decomposition_gsmap.py`,
`15h2_eur_expansion_enrichment.py`, `scripts/figures/fig4a_overview_candidates.py`,
`scripts/figures/fig4{b,e}_snatac_accessibility.R`.

PROVENANCE (recipe locked 2026-07-06 with the PI; see fig4a_overview_candidates.py
docstring). The FINAL set is EXACTLY the union of two arms:

  GENETIC (3,038) — RE-DERIVED HERE from the multi-evidence atlas (S3 genetic
    source). "all coloc (ABF + SuSiE, every GWAS, PP.H4>0.5) + fine-map credible
    sets (SuSiEx / MESuSiE) + regulatory-GWAS-driven (ABC / caQTL / sQTL / cCRE
    credible-set hits) + rare-variant burden + ClinVar P/LP + genome-wide-
    significant TWAS." Encoded as GENETIC_RULES below. This reconstruction
    reproduces the on-disk universe_genetic.txt EXACTLY (3,038/3,038, 0 extra,
    0 missing — validated 2026-07-08).

  TRANSCRIPTOMIC (8,260) — CONSUMED HERE as a locked upstream artifact, NOT
    re-derived. It is the union of "every DE contrast (bulk disease-vs-control,
    strict histologic extremes (definite-disease vs strict-control), MASH-vs-control,
    stage steatosis/SH/cirrhosis, fibrosis-gradient F0→F4, MASH-vs-MASL, sc
    pseudobulk per cell type) gated by the same interval-null TREAT @ lfc=0.25 as
    the core Fig-3 DEG, ∪ hotspot / LIANA / SVG / conserved-core membership."
    (2026-07-10: added the strict-extremes and MASH-vs-control contrasts as an
    explicit dynamic-DEG family; +172 net-new genes over the 2026-07-06 lock of
    8,088 — extremes +163, MASH-vs-control +9. Orthogonality is robust across the
    whole family: 92% canonical / 87% extremes / 97% MASH-vs-ctrl / 98% MASH-vs-MASL,
    see RNA-seq/results/audit_sensitivity/orthogonality_by_contrast.csv.)
    Those per-contrast TREAT refits are owned by the bulk/sc DE pipelines and are
    NOT stored as atlas columns (the atlas keeps padj, not per-contrast treat_fdr),
    so the arm cannot be reproduced from the atlas alone. Source pointers for a
    future full re-derivation:
      - bulk disease-vs-control TREAT core (1,915):
          RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/
          canonical_deg_results.csv  (treat_fdr < 0.05 @ lfc=0.25)
      - stage / fibrosis-gradient / progression contrasts:
          RNA-seq/results/{granular_staging,reversal,stratified_causal}/…
      - MASH-vs-MASL (mash_vs_masl_treat_degs = 207):
          RNA-seq/results/audit_sensitivity/lfc_sweep_mash_vs_masl*.csv
      - strict histologic extremes (definite-disease vs strict-control; TREAT DEGs
          3,548; +163 net-new to tx; added 2026-07-10):
          RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/
          sensitivity/extreme_phenotype_definite_vs_strict.csv  (TREAT reconstructed
          from summary stats; validated vs canonical treat_fdr, Jaccard 1.000)
      - MASH-vs-control (c11; TREAT DEGs 762; +9 net-new to tx; added 2026-07-10):
          RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/
          c11_nash_vs_ctrl_lvqw.csv
      - sc pseudobulk per cell type: Analysis/SingleCell pseudobulk DE outputs
      - membership: hotspot modules, LIANA differential interactions, spatial SVGs,
          conserved-core (spatial_utils.load_conserved)
    The original one-off builder was interactive and was not committed; this file
    documents the recipe so the arm is auditable, and validates the locked list is
    present + non-empty before use.

VALIDATION GATE
---------------
By default the script is READ-ONLY: it recomputes genetic + FINAL, compares to the
on-disk files, and prints a Jaccard report. It writes canonical files ONLY with
--write, and even then refuses to overwrite if the recomputed set differs from the
on-disk lock (writes a `*.rebuilt.txt` sidecar + exits non-zero instead), so drift
can never be laundered into the canonical files silently.

ENV: `spatial` (rnaseq has a numpy/pandas ABI conflict for Python).
USAGE:
    python build_prioritized_universe.py            # validate only (read-only)
    python build_prioritized_universe.py --write    # write canonical files if validated
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS = os.path.join(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
UDIR = os.path.join(BASE, "Analysis/Spatial/results/universe_validation")
F_GEN = os.path.join(UDIR, "universe_genetic.txt")
F_TX = os.path.join(UDIR, "universe_transcriptomic.txt")
F_FINAL = os.path.join(UDIR, "prioritized_universe_FINAL.txt")
F_META = os.path.join(UDIR, "universe_FINAL_meta.json")

# ── GENETIC arm rule set (each rule is a documented subset of the genetic arm) ──
# (column, kind, threshold). kind: 'gt' numeric >thr | 'lt' numeric <thr | 'flag' truthy.
GENETIC_RULES = [
    # COLOC (ABF + SuSiE, every GWAS / ancestry / QTL panel), PP.H4 > 0.5
    ("coloc_abf_best_pp4",            "gt", 0.5),
    ("coloc_best_pp4_polyfun",        "gt", 0.5),
    ("coloc_best_susie_pp4_polyfun",  "gt", 0.5),
    ("coloc_abf_best_pp4_EUR",        "gt", 0.5),
    ("sceqtl_coloc_best_pp4",         "gt", 0.5),
    ("broadaway_coloc_pp4",           "gt", 0.5),
    ("pdff_coloc_pp4",                "gt", 0.5),
    ("best_liver_enzyme_pp4",         "gt", 0.5),
    ("zenodo_nafld_coloc",            "gt", 0.5),
    ("n_coloc_sources",               "gt", 0),   # atlas aggregate: any coloc source (incl ieqtl/zenodo/stage)
    # Fine-mapping credible sets (joint / multi-ancestry)
    ("susiex_n_cs",                   "gt", 0),
    ("mesusie_in_shared_cs",          "flag", None),
    ("mesusie_in_eur_cs",             "flag", None),
    ("mesusie_in_eas_cs",             "flag", None),
    # Regulatory-GWAS-driven (credible-set variant → regulatory element)
    ("abc_masld_gwas_driven",         "flag", None),
    ("caqtl_masld_gwas_driven",       "flag", None),
    ("sqtl_masld_gwas_driven",        "flag", None),
    ("ccre_credset_hit",              "flag", None),
    ("abc_v2g_hit",                   "flag", None),
    ("caqtl_credset_hit",             "flag", None),
    ("sqtl_credset_hit",              "flag", None),
    ("n_regulatory_layers",           "gt", 0),   # atlas aggregate: any regulatory credible-set layer
    # Rare-variant burden + ClinVar + genome-wide-significant TWAS
    ("burden_masld_hit",              "flag", None),
    ("clinvar_n_plp_masld",           "gt", 0),
    ("twas_pval",                     "lt", 2.6e-6),
]


def read_set(path):
    with open(path) as fh:
        return set(x.strip() for x in fh if x.strip())


def apply_rule(atl, sym, col, kind, thr):
    if col not in atl.columns:
        print(f"  [warn] atlas missing column '{col}' — rule skipped", file=sys.stderr)
        return set()
    if kind == "gt":
        return set(sym[pd.to_numeric(atl[col], errors="coerce") > thr])
    if kind == "lt":
        return set(sym[pd.to_numeric(atl[col], errors="coerce") < thr])
    if kind == "flag":
        v = atl[col].astype(str).str.strip().str.lower()
        return set(sym[v.isin(["true", "1", "1.0", "yes"])])
    raise ValueError(kind)


def derive_genetic(atl):
    sym = atl["human_symbol"].astype(str)
    g = set()
    for col, kind, thr in GENETIC_RULES:
        g |= apply_rule(atl, sym, col, kind, thr)
    return g


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a or b) else 1.0


def report(name, got, disk):
    j = jaccard(got, disk)
    extra, miss = got - disk, disk - got
    status = "EXACT" if not extra and not miss else "DRIFT"
    print(f"  {name:14s} got={len(got):6d} disk={len(disk):6d} "
          f"Jaccard={j:.4f} extra={len(extra)} missing={len(miss)}  [{status}]")
    if extra:
        print(f"      extra (got−disk) sample: {sorted(extra)[:12]}")
    if miss:
        print(f"      missing (disk−got) sample: {sorted(miss)[:12]}")
    return not extra and not miss


def main():
    ap = argparse.ArgumentParser(description="Rebuild + validate the prioritized universe.")
    ap.add_argument("--write", action="store_true",
                    help="write canonical files (only if recomputed sets match the on-disk lock)")
    args = ap.parse_args()

    for p in (ATLAS, F_GEN, F_TX, F_FINAL):
        if not os.path.exists(p):
            sys.exit(f"[fatal] required input missing: {p}")

    print(f"[build-universe] atlas: {ATLAS}")
    atl = pd.read_csv(ATLAS, low_memory=False)

    # GENETIC — re-derived from atlas
    gen_new = derive_genetic(atl)
    gen_disk = read_set(F_GEN)

    # TRANSCRIPTOMIC — locked upstream artifact (validate, do not re-derive)
    tx_disk = read_set(F_TX)
    if not tx_disk:
        sys.exit(f"[fatal] transcriptomic lock is empty: {F_TX}")

    # FINAL — union
    final_new = gen_new | tx_disk
    final_disk = read_set(F_FINAL)

    print("\n[validation] recomputed vs on-disk lock")
    ok_gen = report("genetic", gen_new, gen_disk)
    ok_final = report("FINAL", final_new, final_disk)
    # sanity: on-disk FINAL must itself equal disk(gen) ∪ disk(tx)
    ok_union = (gen_disk | tx_disk) == final_disk
    print(f"  {'union-check':14s} disk(gen)∪disk(tx) == disk(FINAL): {ok_union}")
    print(f"  transcriptomic: locked upstream artifact, {len(tx_disk)} genes "
          f"(NOT re-derived — see header for source pointers)")

    all_ok = ok_gen and ok_final and ok_union
    print(f"\n[validation] overall: {'PASS' if all_ok else 'FAIL'}")

    meta = {
        "universe_final": len(final_new),
        "genetic": len(gen_new),
        "transcriptomic": len(tx_disk),
        "genetic_rederived_exact": ok_gen,
        "final_validated_exact": ok_final,
        "genetic_rules": [f"{c} {k} {t}" for c, k, t in GENETIC_RULES],
        "transcriptomic_source": "locked upstream artifact (per-contrast TREAT + hotspot/LIANA/SVG/conserved-core)",
        "recipe_locked": "2026-07-06",
        "recipe_updated": "2026-07-10 (+strict-extremes, +MASH-vs-control; tx 8,088->8,260, FINAL 9,882->10,044)",
        "tx_augment_2026_07_10": {"strict_extremes_treat_degs": 3548,
                                  "mash_vs_control_treat_degs": 762,
                                  "net_new_to_tx": 172},
        "rebuilt_by": "Analysis/Spatial/scripts/build_prioritized_universe.py",
    }
    # carry forward documented tx-arm provenance counts if the existing meta has them
    if os.path.exists(F_META):
        try:
            old = json.load(open(F_META))
            for k in ("mash_vs_masl_treat_degs", "mvm_new", "validated_orthogonal", "proteo", "scatac"):
                if k in old:
                    meta[k] = old[k]
        except Exception as e:
            print(f"  [warn] could not read existing meta: {e}", file=sys.stderr)

    if not args.write:
        print("\n[dry-run] read-only. Re-run with --write to update canonical files "
              "(gated on validation PASS).")
        return 0 if all_ok else 1

    # --write path: refuse to clobber on drift
    if not (ok_gen and ok_final):
        side_g = F_GEN.replace(".txt", ".rebuilt.txt")
        side_f = F_FINAL.replace(".txt", ".rebuilt.txt")
        with open(side_g, "w") as fh:
            fh.write("\n".join(sorted(gen_new)) + "\n")
        with open(side_f, "w") as fh:
            fh.write("\n".join(sorted(final_new)) + "\n")
        sys.exit(f"[refused] recomputed sets differ from the on-disk lock — wrote sidecars "
                 f"instead of overwriting:\n  {side_g}\n  {side_f}\n"
                 f"Reconcile before promoting to canonical.")

    with open(F_GEN, "w") as fh:
        fh.write("\n".join(sorted(gen_new)) + "\n")
    with open(F_FINAL, "w") as fh:
        fh.write("\n".join(sorted(final_new)) + "\n")
    with open(F_META, "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"\n[write] updated (validated exact):\n  {F_GEN}\n  {F_FINAL}\n  {F_META}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
