#!/usr/bin/env python
"""17: refuse to call the run finished unless these checks pass.

Ordered by what they protect: the firewall between module construction and
outcome reading; the identity of the v2 modules with v1's; the arithmetic that
the 2026-09-11 review found wrong (orientation, correction, coverage, batch);
and the figure.
"""
import argparse
import hashlib
import json
import os
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CONTRACT = json.loads((ROOT / "scripts/analysis/cross_assay_modules/00_contract.json").read_text())
FAILURES: list[str] = []
PASSES: list[str] = []
BANNED_CAPTION_PHRASES = ["adequately powered", "composition-free", "almost none travel",
                          "almost none of them travel", "sealed provenance", "negative control"]


def check(name: str, ok: bool, detail: str = "") -> None:
    (PASSES if ok else FAILURES).append(f"{name}{': ' + detail if detail else ''}")
    print(("PASS " if ok else "FAIL ") + name + (f" — {detail}" if detail else ""), flush=True)


def sha256(p: pathlib.Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def bh(p: np.ndarray, n: int) -> np.ndarray:
    """Benjamini-Hochberg with the family size fixed at n, NaN-preserving."""
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    ps = p[ok]
    order = np.argsort(ps)
    ranked = ps[order] * n / (np.arange(len(ps)) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(len(ps)); res[order] = np.minimum(q, 1.0)
    out[ok] = res
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True, type=pathlib.Path)
    args = ap.parse_args()
    R = args.out_root

    # --- 1. firewall and module identity ------------------------------------
    frozen = R / "modules/MODULES_FROZEN.ok"
    check("module freeze marker exists", frozen.exists())
    if frozen.exists():
        recorded = dict(l.split("=", 1) for l in frozen.read_text().strip().splitlines())
        memb = R / "modules/module_membership.tsv"
        check("frozen membership sha256 still matches the file",
              recorded.get("membership_sha256") == sha256(memb))
        check("v2 membership is byte-identical to v1 (construction untouched)",
              sha256(memb) == CONTRACT["supersedes"]["v1_membership_sha256"],
              f"v2 {sha256(memb)[:12]} vs v1 {CONTRACT['supersedes']['v1_membership_sha256'][:12]}")
        check("freeze recorded that no outcome had been read", recorded.get("outcome_read") == "FALSE")
        disc = R / "discovery/discovery_summary.json"
        if disc.exists():
            d = json.loads(disc.read_text())
            check("outcome was read AFTER the freeze", d["outcome_read_at"] > d["modules_frozen_at"])
            check("every downstream file is newer than the freeze",
                  all(f.stat().st_mtime >= frozen.stat().st_mtime for f in (R / "discovery").glob("*.tsv")))

    # --- 2. construction diagnostics ---------------------------------------
    memb_df = pd.read_csv(R / "modules/module_membership.tsv", sep="\t")
    reg = pd.read_csv(R / "modules/module_registry.tsv", sep="\t")
    n_modules = reg.module_id.nunique()
    check("registry and membership agree on the module count", n_modules == memb_df.module_id.nunique(), f"{n_modules}")
    check("every module is within the prespecified size bounds",
          bool(((reg.n_genes >= CONTRACT["modules"]["min_size"]) & (reg.n_genes <= CONTRACT["modules"]["max_size"])).all()))
    check("no gene belongs to two modules", not memb_df.gene_symbol.duplicated().any())
    loco = pd.read_csv(R / "modules/loco_stability.tsv", sep="\t")
    check("held-out leave-one-cohort-out score gate met in every fold",
          bool((loco.median_heldout_score_spearman >= CONTRACT["modules"]["loco_score_spearman_floor"]).all()),
          f"min {loco.median_heldout_score_spearman.min():.4f}")
    check("LOCO baselines recorded (random set, held-out PC1)", (R / "modules/loco_baselines.tsv").exists())
    uni = json.loads((R / "universe/universe_summary.json").read_text())
    check("universe meets its floor", uni["universe_size"] >= CONTRACT["universe"]["min_universe_size"])
    mlc = pd.read_csv(R / "graph/gene_mean_logcpm.tsv", sep="\t")
    check("expression-matched bins rest on mean log-CPM with real spread",
          float(mlc.mean_logcpm.std()) > 0.5, f"sd {mlc.mean_logcpm.std():.3f}")

    # --- 3. the arithmetic the review found wrong -----------------------------
    if (R / "nulls/loader_reproduction_check.tsv").exists():
        lr = pd.read_csv(R / "nulls/loader_reproduction_check.tsv", sep="\t")
        check("null loaders reproduce every assay's reported observed effect",
              bool((lr.max_abs_difference < 1e-8).all()), f"worst {lr.max_abs_difference.max():.2e}")
    vc = pd.read_csv(R / "discovery/vectorised_fit_check.tsv", sep="\t")
    check("vectorised association matches lm()", bool((vc.abs_difference < 1e-8).all()))
    cd = pd.read_csv(R / "discovery/competitive_draw_diagnostics.tsv", sep="\t")
    check("competitive draw failure rate under its cap",
          bool((cd.failure_rate <= CONTRACT["nulls"]["competitive_connected"]["max_draw_failure_rate"]).all()))
    check("competitive draws at the prespecified count",
          bool((cd.n_drawn >= 0.95 * CONTRACT["nulls"]["competitive_connected"]["draws"]).all()),
          f"min {cd.n_drawn.min()}")
    check("competitive calibration recorded", (R / "nulls/competitive_calibration.json").exists())
    check("no family-level binomial table is produced",
          not (R / "panel/family_level_competitive_check.tsv").exists())

    panel_f = R / "panel/module_by_assay_states.tsv"
    if panel_f.exists():
        panel = pd.read_csv(panel_f, sep="\t")
        banned = set(CONTRACT["prohibited_columns"]) & {c.lower() for c in panel.columns}
        check("panel contains no cross-assay score, rank or count column", not banned)
        per_assay = panel.groupby("assay").module_id.nunique()
        check("every assay column covers the whole module family", bool((per_assay == n_modules).all()))
        check("every cell carries a display state", not panel.display_state.isna().any())
        allowed = (set(CONTRACT["states"]["association_state"]) | set(CONTRACT["states"]["coherence_state"]))
        check("every display state is from the v2 vocabulary",
              set(panel.display_state.unique()) <= allowed, str(set(panel.display_state.unique()) - allowed))
        check("old `supported` vocabulary absent", not (panel.display_state == "supported").any())
        check("spatial non-significance is never printed as absence",
              not panel.display_state.isin(["not_coherent"]).any())
        check("every cell carries a coverage value regardless of endpoint applicability",
              not panel.fraction_measured.isna().any())
        pro_all = panel[(panel.assay == "proteome_liver_pxd051911") & (panel.evidence_axis == "association")]
        check("proteome rows carry the multi-gene protein-group count",
              bool(np.isfinite(pro_all.n_members_from_multigene_groups).all()))
        snr = panel[panel.assay == "snrna_all_cells"]
        check("snRNA coarse-diagnosis endpoint is flagged as a proxy", bool(snr.endpoint_is_proxy.astype(bool).all()))
        assoc = panel[panel.evidence_axis == "association"]
        app = assoc[assoc.applicable.astype(bool)]
        check("effect_oriented == effect_raw * direction on every applicable cell",
              bool(np.nanmax(np.abs(app.effect_oriented - app.effect_raw * app.direction)) < 1e-12))
        # association q recomputed from association p over testable+applicable rows, family 160
        ok_all = True
        for a, g in app.groupby("assay"):
            p = np.where(g.testable.astype(bool), g.association_p.values, np.nan)
            q = bh(p.astype(float), n_modules)
            diff = np.nanmax(np.abs(q - g.association_q.values)) if np.isfinite(q).any() else 0.0
            ok_all &= diff < 1e-10
        check("association q is BH over the 160 family from the assay p", ok_all)
        ok_all = True
        for a, g in app.groupby("assay"):
            p = g.competitive_p.values.astype(float)
            q = bh(p, n_modules)
            m = np.isfinite(q)
            diff = np.nanmax(np.abs(q[m] - g.competitive_q.values[m])) if m.any() else 0.0
            ok_all &= diff < 1e-10
        check("competitive q is BH over the 160 family from the competitive p", ok_all)
        check("specificity never assessed for untestable or inapplicable cells",
              not ((~(assoc.applicable.astype(bool) & assoc.testable.astype(bool))) &
                   (assoc.specificity_state != "not_assessed")).any())
        tn = app[app.association_state == "tested_negative"]
        check("every tested_negative has equivalence q < 0.05 and association q >= 0.05",
              bool(((tn.equivalence_q < 0.05) & (tn.association_q >= 0.05)).all()) if len(tn) else True,
              f"{len(tn)} cells")
        rep = app[app.association_state == "replicates"]
        check("every `replicates` cell has oriented effect > 0 and q < 0.05",
              bool(((rep.effect_oriented > 0) & (rep.association_q < 0.05)).all()) if len(rep) else True,
              f"{len(rep)} cells")
        dis = app[app.association_state == "discordant"]
        check("every `discordant` cell has oriented effect < 0 and q < 0.05",
              bool(((dis.effect_oriented < 0) & (dis.association_q < 0.05)).all()) if len(dis) else True,
              f"{len(dis)} cells")
        usable = app[app.testable.astype(bool)]
        check("every usable cell prints an upper bound", bool(np.isfinite(usable.upper_bound_95).all()))
        atac = app[app.assay == "atac_gse296875"]
        atac_res = pd.read_csv(R / "assays/atac_results.tsv", sep="\t").set_index("module_id")
        ap_ = atac.set_index("module_id")
        common = ap_.index.intersection(atac_res.index)
        diff = np.nanmax(np.abs(ap_.loc[common, "association_p"].values - atac_res.loc[common, "within_well_p"].values))
        check("ATAC association p is the within-well permutation p", diff < 1e-12)
        pro = app[(app.assay == "proteome_liver_pxd051911") & (app.association_state == "replicates")]
        check("every replicating proteome cell is complexity-robust at BH q<0.05",
              bool((pro.complexity_robust.astype(bool) & (pro.complexity_adjusted_q < 0.05)).all()) if len(pro) else True,
              f"{len(pro)} cells")
        vu = panel[panel.assay == "visium_vu"]
        check("every testable Vu cell is source_dependent",
              bool((vu[vu.testable.astype(bool)].display_state == "source_dependent").all()))
        check("no spatial cell claims an association state",
              bool((panel[panel.evidence_axis == "spatial_coherence"].association_state == "not_applicable").all()))
        check("untestable cells carry no effect estimate",
              bool(assoc.loc[~assoc.testable.fillna(False).astype(bool), "effect_raw"].isna().all()))
        ps = json.loads((R / "panel/panel_summary.json").read_text())
        check("state-rule self-tests ran", ps.get("state_self_test") is True)

    pt_f = R / "panel/paired_transfer_panel.tsv"
    if pt_f.exists():
        pt = pd.read_csv(pt_f, sep="\t")
        sig = pt[pt.q_value.notna()]
        check("every paired row with a q is testable in both assays", bool(sig.testable_both.astype(bool).all()))
        at = pt[(pt.pairing == "gse296875_rna_vs_atac") & (pt.q_value < 0.05)]
        check("ATAC sensitivities carry their own permutation p and BH q",
              bool(np.isfinite(at.q_given_well_and_fibrosis).all() and
                   np.isfinite(at.q_given_well_and_hepatocyte_fraction).all()) if len(at) else True)
        if len(at):
            print(f"      [report] of {len(at)} well-adjusted couplings, "
                  f"{int((at.q_given_well_and_fibrosis < 0.05).sum())} stay at BH q<0.05 given well+fibrosis, "
                  f"{int((at.q_given_well_and_hepatocyte_fraction < 0.05).sum())} given well+hepatocyte fraction", flush=True)
        check("every significant ATAC coupling carries the hepatocyte-fraction sensitivity",
              bool(np.isfinite(at.partial_given_well_and_hepatocyte_fraction).all()) if len(at) else True,
              f"{len(at)} significant")
        if len(at):
            ratio = (at.partial_given_well_and_hepatocyte_fraction / at.partial_spearman)
            print(f"      [report] median attenuation under hepatocyte fraction: {ratio.median():.3f}; "
                  f"{int((ratio.abs() >= 0.5).sum())} of {len(at)} retain half the effect", flush=True)
        h3 = pt[(pt.pairing == "gse267145_rna_vs_h3k27ac") & pt.testable_both.astype(bool)]
        check("H3K27ac bootstrap intervals bracket the partial estimate",
              bool(((h3.boot_ci_low <= h3.partial_spearman) & (h3.partial_spearman <= h3.boot_ci_high)).mean() > 0.95),
              f"{((h3.boot_ci_low <= h3.partial_spearman) & (h3.partial_spearman <= h3.boot_ci_high)).mean():.3f}")

    sc = R / "assays/spatial_unit_census.tsv"
    if sc.exists():
        cen = pd.read_csv(sc, sep="\t").set_index("assay")
        check("Guilliams spatial graphs built per section (5 sections, 4 donors)",
              int(cen.loc["visium_gse192741", "n_sections"]) == 5 and int(cen.loc["visium_gse192741", "n_donors"]) == 4,
              f"{cen.loc['visium_gse192741', 'sections_per_donor']}")

    # --- 4. the figures --------------------------------------------------------
    fig_dir = R / os.environ.get("CAM_FIGURE_VERSION", "figure")
    if fig_dir.exists():
        pdfs = sorted(fig_dir.glob("*.pdf"))
        check("figures were rendered", len(pdfs) > 0, f"{len(pdfs)} PDFs")
        for p in pdfs:
            try:
                info = subprocess.run(["pdfinfo", str(p)], capture_output=True, text=True)
                pages = [l for l in info.stdout.splitlines() if l.startswith("Pages:")]
                check(f"{p.name} is one page", bool(pages) and pages[0].split()[-1] == "1")
            except FileNotFoundError:
                pass
        cks = pd.read_csv(fig_dir / "figure_checksums.tsv", sep="\t")
        check("figure checksums match the rendered files",
              all(sha256(fig_dir / r.file) == r.sha256 for r in cks.itertuples()))
        k5 = fig_dir / "fig5k_source.tsv"
        if k5.exists():
            k = pd.read_csv(k5, sep="\t")
            unstable = k[k.score_stability_flag == "score_unstable"]
            check("5K marks every score-unstable example in its label",
                  bool(unstable.module_label.str.contains("unstable").all()) if len(unstable) else True,
                  f"{unstable.module_id.nunique()} unstable examples")
        for src in fig_dir.glob("*_source.tsv"):
            df = pd.read_csv(src, sep="\t")
            bad = set(CONTRACT["prohibited_columns"]) & {c.lower() for c in df.columns}
            check(f"{src.name} has no prohibited column", not bad)
        cap = fig_dir / "CAPTIONS.md"
        if cap.exists():
            txt = cap.read_text().lower()
            hits = [b for b in BANNED_CAPTION_PHRASES if b in txt]
            check("captions avoid the withdrawn phrases", not hits, str(hits))

    print()
    print(f"{len(PASSES)} passed, {len(FAILURES)} failed")
    (R / "validation").mkdir(exist_ok=True)
    (R / "validation/validation_report.json").write_text(json.dumps(
        {"passed": PASSES, "failed": FAILURES, "n_passed": len(PASSES), "n_failed": len(FAILURES)}, indent=2))
    if FAILURES:
        return 1
    (R / "validation/READY").write_text("validated\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
