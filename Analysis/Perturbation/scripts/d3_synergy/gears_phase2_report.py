#!/usr/bin/env python
"""Phase E — compare Phase 1 (Norman K562 zs) vs Phase 2 (Saunders ft) GEARS
on the same D3 invocab + Saunders-internal pair sets.

Outputs:
  Analysis/Perturbation/results/d3_synergy/gears_v2_finetune_report.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
D3_RESULTS = PROJECT_ROOT / "Analysis/Perturbation/results/d3_synergy"


def load_envelope(path):
    if not path.exists():
        return None
    with open(path) as fh:
        return json.load(fh)


def summarize(env):
    if env is None:
        return dict(n=0, classes={}, sigma_mean=np.nan, sigma_max=np.nan, k562_bias=np.nan)
    preds = env.get("predictions", [])
    if not preds:
        return dict(n=0, classes={}, sigma_mean=np.nan, sigma_max=np.nan,
                    k562_bias=env.get("notes", ""))
    df = pd.DataFrame(preds)
    return dict(
        n=len(df),
        classes=df["synergy_class"].value_counts().to_dict() if "synergy_class" in df else {},
        sigma_mean=float(df["sigma_above_additive"].mean()) if "sigma_above_additive" in df else np.nan,
        sigma_max=float(df["sigma_above_additive"].max()) if "sigma_above_additive" in df else np.nan,
        k562_bias=float(df["k562_bias_confidence"].mean()) if "k562_bias_confidence" in df else np.nan,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-md", default=str(D3_RESULTS / "gears_v2_finetune_report.md"))
    ap.add_argument("--ckpt-dir", default="data/perturbation/finetuned/gears/saunders_hep_v1")
    args = ap.parse_args()

    phase1 = load_envelope(D3_RESULTS / "gears_zero_shot_tier1_pairs__all.json")
    p2_tier1 = load_envelope(D3_RESULTS / "gears_v2_finetune_tier1_pairs__invocab.json")
    p2_tier2 = load_envelope(D3_RESULTS / "gears_v2_finetune_tier2_pairs__invocab.json")
    p2_saund = load_envelope(D3_RESULTS / "gears_v2_finetune_saunders_internal__top200.json")

    train_meta_path = PROJECT_ROOT / args.ckpt_dir / "train_metadata.json"
    train_meta = json.loads(train_meta_path.read_text()) if train_meta_path.exists() else {}

    lines = []
    lines.append("# GEARS Phase 2 Fine-tune Report")
    lines.append("")
    lines.append("**Substrate**: Saunders 2025 Perturb-Multi mouse hepatocyte "
                 "(GSE275673; 203 mouse target genes / 181 1:1 human orthologs / "
                 "~103K confident-tier cells).")
    lines.append("**Base**: Norman 2019 K562 essential (combinatorial Perturb-seq) — same Phase 1 substrate; no warm-start in this run (vocab mismatch precludes weight transfer without union-vocab engineering).")
    lines.append("")

    lines.append("## Training summary")
    lines.append("")
    if train_meta:
        a = train_meta.get("args", {})
        lines.append(f"- cells trained: {train_meta.get('n_cells', '?'):,}")
        lines.append(f"- gene-expression universe: {train_meta.get('n_genes', '?'):,}")
        lines.append(f"- conditions (incl ctrl): {train_meta.get('n_conditions', '?')}")
        lines.append(f"- perturbable vocab (`pert_names`): {train_meta.get('pert_universe_size', '?')}")
        lines.append(f"- control cells: {train_meta.get('control_cells', '?'):,}")
        lines.append(f"- epochs={a.get('epochs','?')}, lr={a.get('lr','?')}, "
                     f"batch_size={a.get('batch_size','?')}, hidden_size={a.get('hidden_size','?')}")
        lines.append(f"- holdout_frac={a.get('holdout_frac','?')}")
        lines.append(f"- training device: {train_meta.get('device','?')}")
        lines.append(f"- training time: {train_meta.get('training_time_sec', 0)/60:.1f} min")
    else:
        lines.append("(train_metadata.json not found — fine-tune may not have completed)")
    lines.append("")

    lines.append("## Phase 1 (Norman K562 zero-shot) vs Phase 2 (Saunders fine-tuned)")
    lines.append("")
    lines.append("| metric | Phase 1 (Norman K562 zs) | Phase 2 tier1 | Phase 2 tier2 | Phase 2 Saunders-internal top200 |")
    lines.append("|---|---|---|---|---|")
    p1s = summarize(phase1)
    p2t1s = summarize(p2_tier1)
    p2t2s = summarize(p2_tier2)
    p2sas = summarize(p2_saund)
    lines.append(f"| n_predictions | {p1s['n']} | {p2t1s['n']} | {p2t2s['n']} | {p2sas['n']} |")
    lines.append(f"| sigma_mean | {p1s['sigma_mean']:.3f} | {p2t1s['sigma_mean']:.3f} | "
                 f"{p2t2s['sigma_mean']:.3f} | {p2sas['sigma_mean']:.3f} |")
    lines.append(f"| sigma_max | {p1s['sigma_max']:.3f} | {p2t1s['sigma_max']:.3f} | "
                 f"{p2t2s['sigma_max']:.3f} | {p2sas['sigma_max']:.3f} |")
    lines.append(f"| k562_bias_confidence | {p1s['k562_bias']:.2f} | "
                 f"{p2t1s['k562_bias'] if isinstance(p2t1s['k562_bias'], float) else 'n/a'} | "
                 f"{p2t2s['k562_bias'] if isinstance(p2t2s['k562_bias'], float) else 'n/a'} | "
                 f"{p2sas['k562_bias'] if isinstance(p2sas['k562_bias'], float) else 'n/a'} |")
    lines.append("")

    lines.append("## Key findings")
    lines.append("")
    lines.append("- **Vocab constraint**: After Saunders fine-tune, GEARS' perturbable universe = 181 genes "
                 "(the human orthologs of Saunders' 203 mouse target genes; 1:1 orthologs only).")
    lines.append("- **Tier1/Tier2 invocab pairs**: 0 / 4 tier1 unique genes (TADA1, VPS41, HECTD1, RNF103) "
                 "are in Saunders' panel; 6 / 81 tier2 unique genes are (GPN1, MTOR, RPL8, SLU7, SNRPA1, SRP72) "
                 "but no tier2 PAIR has both genes in Saunders. So Phase 2 cannot improve on the Phase 1 "
                 "null result for these specific invocab hits — both Phase 1 (Norman K562 vocab miss) and "
                 "Phase 2 (Saunders vocab miss) return 0 predictions for these targets.")
    lines.append("- **Saunders-internal pair set**: Among C(181,2)=16,290 in-vocab pairs, the top 200 ranked "
                 "by `|bulk_tstat| + 5×SuSiE-COLOC PP4` capture liver-relevant synergy candidates "
                 "(MET × MTOR, INSR × MET, FASN × MTOR, …). These are predictable by the Saunders-tuned "
                 "model and are the proper Phase 2 deliverable given the vocab constraint.")
    lines.append("")

    if p2sas['n'] > 0 and p2_saund is not None:
        df = pd.DataFrame(p2_saund["predictions"])
        df_sorted = df.sort_values("synergy_magnitude", ascending=False).head(10)
        lines.append("## Top-10 synergistic Saunders-internal pairs (Phase 2 fine-tuned)")
        lines.append("")
        lines.append("| gene1 | gene2 | sigma_above_additive | synergy_magnitude | synergy_class |")
        lines.append("|---|---|---|---|---|")
        for _, r in df_sorted.iterrows():
            lines.append(f"| {r['gene1']} | {r['gene2']} | "
                         f"{r['sigma_above_additive']:.3f} | "
                         f"{r['synergy_magnitude']:.3f} | {r['synergy_class']} |")
        lines.append("")

    lines.append("## K562-bias-confidence shift")
    lines.append("")
    lines.append(f"- Phase 1 (Norman K562 zs): hard-coded `k562_bias_confidence = 1.0` "
                 "(K562 lineage, no hepatocyte signal).")
    lines.append(f"- Phase 2 (Saunders ft): hard-coded `k562_bias_confidence = 0.2` "
                 "(mouse hepatocyte lineage; minimal K562 transfer since no warm-start). "
                 "If union-vocab fine-tune (Norman ∪ Saunders) is requested in a follow-up, "
                 "k562_bias would land closer to 0.5.")
    lines.append("")

    lines.append("## Files produced")
    lines.append("")
    lines.append(f"- Fine-tuned checkpoint: `{args.ckpt_dir}/best.ckpt/`")
    lines.append(f"- Prepped AnnData: `{args.ckpt_dir}/prepped.h5ad`")
    lines.append(f"- PertData processed dataset: `{args.ckpt_dir}/saunders_hep/`")
    lines.append("- Predictions:")
    lines.append(f"  - `Analysis/Perturbation/results/d3_synergy/gears_v2_finetune_tier1_pairs__invocab.json` (tier1)")
    lines.append(f"  - `Analysis/Perturbation/results/d3_synergy/gears_v2_finetune_tier2_pairs__invocab.json` (tier2)")
    lines.append(f"  - `Analysis/Perturbation/results/d3_synergy/gears_v2_finetune_saunders_internal__top200.json` (Saunders-internal)")
    lines.append(f"- Hit lists for Saunders-vocab targets:")
    lines.append(f"  - `Analysis/Perturbation/data/hits/d3_saunders_internal_pairs.csv` (16,290 pairs)")
    lines.append(f"  - `Analysis/Perturbation/data/hits/d3_saunders_internal_pairs_top200.csv` (top 200)")
    lines.append(f"  - `Analysis/Perturbation/data/hits/d3_saunders_masld_curated_pairs.csv` (820 MASLD-curated)")

    out = Path(args.out_md)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
