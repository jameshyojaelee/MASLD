#!/usr/bin/env python3
"""Stage D 2026-11-18 captions: reads exclusively from out_20261118/stageD/tables/."""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

D = L.PKG / "stageD"
OUT = L.PKG / "out_20261118/stageD"
TAB = OUT / "tables"
PAN = OUT / "panels"
COND = {
    "bulk_liver": "bulk liver (UBERON:0002107)",
    "hepatocyte": "hepatocyte (UBERON:0002107_CL:0000182)",
    "monocyte_ctrl": "CD14 monocyte control (UBERON:0000178_CL:0001054)",
}


def f3(x):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{float(x):+.3f}"


def txt(x):
    return "" if (x is None or (isinstance(x, float) and np.isnan(x))) else str(x)


def safe(u):
    return u.replace(":", "_")


def main():
    v = pd.read_csv(D / "inputs/stageD_variants.tsv", sep="\t")
    me = pd.read_csv(D / "measured/measured_at_variants.tsv", sep="\t").set_index("uid")
    mc = pd.read_csv(D / "motifs/motif_changes_by_cluster.tsv", sep="\t")
    loc = pd.read_csv(TAB / "stageD_local_change.tsv", sep="\t")
    ism = pd.read_csv(TAB / "stageD_ism_summary.tsv", sep="\t")
    mot_path = TAB / "stageD_motif_ism_overlap.tsv"
    mot = pd.read_csv(mot_path, sep="\t") if mot_path.stat().st_size > 1 else pd.DataFrame(columns=["uid"])
    sg = pd.read_csv(TAB / "stageD_measured_sign_agreement.tsv", sep="\t").set_index("uid")
    for _, r in v.iterrows():
        u = r.uid
        gene = str(r.genes_stageB).split(",")[0].strip() if isinstance(r.genes_stageB, str) else ""
        folder = PAN / f"{r.locus_id}_{gene}".replace(" ", "")
        folder.mkdir(parents=True, exist_ok=True)
        out = []
        out.append(f"# {u} (GRCh38; hg19 {r.hg19_key})\n")
        out.append(
            "Predicted by AlphaGenome 2 (ModelVersion.ALL_FOLDS; candidate illustration, not adopted, not evidence of a causal variant or gene).\n"
        )
        if r.role == "control_SORT1":
            out.append(
                "Role: positive-control variant at the SORT1/CELSR2/PSRC1 locus (rs12740374; Musunuru et al. 2010, Nature 466:714, PMID 20686566).\n"
            )
        else:
            out.append(
                f"Role: selected tier-1 member. Max PIP {float(r.max_pip):.3f} in {r.max_pip_set} "
                f"({r.max_pip_study}, {r.max_pip_ancestry}, {r.max_pip_method}); {int(r.n_sets)} credible sets; "
                f"consequence {r.consequence}. Stage B sign: {r.stageB_sign}.\n"
            )
        out.append("## Panels\n")
        out.append(
            f"- `{safe(u)}_tracks_<condition>_<output>.pdf`: observed-track mean predicted signal, +/-1 kb around the variant "
            "(1-bp resolution); `_zoom150` files show +/-150 bp."
        )
        out.append(
            f"- `{safe(u)}_ism_ref_<condition>_<output>.pdf`: in silico mutagenesis on the REF sequence, 64 bp centred on the variant."
        )
        out.append(f"- `{safe(u)}_ism_alt_hepatocyte_SC_ATAC.pdf`: the same on the ALT-allele background, 32 bp.\n")
        out.append("## Predicted local change (ALT minus REF, log2(sum + 1) over ChromatinScorer window: 501 bp ATAC/DNASE, 1201 bp SC_ATAC)\n")
        out.append("| condition | output | window_bp | observed tracks | Stage D local change | Stage B raw | Stage B percentile |")
        out.append("|---|---|---|---|---|---|---|")
        for _, x in loc[loc.uid == u].iterrows():
            out.append(
                f"| {COND.get(x.condition, x.condition)} | {x.output} | {int(x.window_bp)} | {int(x.n_observed_tracks)} | "
                f"{f3(x.local_log2_alt_minus_ref_matched)} | {f3(x.stageB_raw)} | "
                f"{'NA' if not np.isfinite(x.stageB_percentile) else f'{x.stageB_percentile:.3f}'} |"
            )
        out.append("\n## Mutagenesis summary\n")
        out.append("| ISM | condition | output | variant-position rank by abs importance | positions | own-substitution score | max abs substitution |")
        out.append("|---|---|---|---|---|---|---|")
        for _, x in ism[ism.uid == u].iterrows():
            out.append(
                f"| {x.step} | {x.condition} | {x.output} | {int(x.variant_position_rank_abs_importance)} | {int(x.of_positions)} | "
                f"{f3(x.own_substitution_score)} | {x.max_abs_substitution:.3f} |"
            )
        out.append("\nTop-importance runs (position relative to variant; bases; mean importance; A3 motif hits overlapping):\n")
        for _, x in ism[ism.uid == u].iterrows():
            out.append(f"- {x.step} {x.condition} {x.output}: {x.top_positions}")
        out.append("\n## Motif changes at the variant\n")
        mm = mc[mc.uid == u]
        if not len(mm):
            out.append("No motif cluster gained or lost.")
        for _, m in mm.iterrows():
            out.append(
                f"- {m['name']} ({m.cluster_family}): {m.direction}, delta {m.delta_alt_minus_ref:+.2f} "
                f"(REF {m.s_ref:.2f}, ALT {m.s_alt:.2f}, threshold {m.threshold:.3f}); "
                f"TF expressed in hepatocytes: {m.tf_expressed_hepatocyte}."
            )
        mo = mot[mot.uid == u] if len(mot) else mot
        for _, m in mo.iterrows():
            out.append(
                f"  - {m['name']} {m.direction}: hepatocyte {m.ism} mean abs importance inside hit {m.hit_span_rel} "
                f"vs outside = {('NA' if not np.isfinite(m.importance_ratio_inside_vs_outside) else f'{m.importance_ratio_inside_vs_outside:.2f}')}x."
            )
        out.append("\n## Measured data (not model output)\n")
        e = me.loc[u] if u in me.index else None
        s = sg.loc[u]
        ca = txt(s.caqtl_beta_alt) if e is not None else ""
        out.append(f"- Currin 2025 liver caQTL beta (hg38 ALT): {ca or 'not measured'}; sign agreement with predicted bulk-liver ATAC: {s.caqtl_sign_agrees}.")
        if e is not None and np.isfinite(e.get("hep_ai_pooled_log2_alt_ref", np.nan)):
            out.append(
                f"- Hepatocyte allelic imbalance (D1, pooled): log2 ALT/REF {e.hep_ai_pooled_log2_alt_ref:+.3f}, "
                f"z {e.hep_ai_pooled_z:+.2f}, {int(e.hep_ai_pooled_n_donors)} donors; sign agreement with predicted hepatocyte SC_ATAC: {s.hep_ai_sign_agrees}."
            )
        else:
            out.append("- Hepatocyte allelic imbalance: not measured.")
        out.append(f"- MPRA: {txt(e.mpra) if e is not None else 'not tested'}.")
        sus = txt(e.coloc_r2_susie_genes) if e is not None else ""
        abf = txt(e.coloc_r2_abf_genes) if e is not None else ""
        out.append(f"- Colocalized eQTL genes at the locus (COLOC r2): SuSiE {sus or 'none'}; ABF {abf or 'none'}.")
        (folder / f"CAPTION_{safe(u)}.md").write_text("\n".join(out) + "\n")
    L.log(f"Stage D 20261118 captions written for {len(v)} variants into {PAN}")


if __name__ == "__main__":
    main()
