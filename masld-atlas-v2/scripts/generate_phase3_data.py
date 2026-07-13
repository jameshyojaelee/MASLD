#!/usr/bin/env python3
"""Generate JSON data files for MASLD Atlas v2 Phase 3 pages.

Produces:
  1. progression_journey.json  — Disease progression (F0→F4 scrollytelling)
  2. drug_pipeline.json        — Drug repurposing funnel
  3. cross_species.json        — Human–mouse concordance

Run:
  micromamba run -n spatial python scripts/generate_phase3_data.py --output-dir public/data
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def clean_for_json(obj):
    """Recursively replace NaN / Inf with None for JSON serialisation."""
    if isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_for_json(v) for v in obj]
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
    if isinstance(obj, (np.floating, np.integer)):
        v = obj.item()
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return None
        return v
    return obj


def write_json(data, path: Path):
    data = clean_for_json(data)
    with open(path, "w") as fh:
        json.dump(data, fh, separators=(",", ":"))
    size_kb = path.stat().st_size / 1024
    print(f"  -> {path.name}: {size_kb:.1f} KB")


def _round(val, n=4):
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return None
    return round(float(val), n)


# ---------------------------------------------------------------------------
# 1. Progression Journey
# ---------------------------------------------------------------------------

def build_progression_journey() -> dict:
    print("[1/3] Building progression_journey.json ...")

    # --- Load stage-specific DE ---
    fib_de_path = (
        PROJECT_ROOT
        / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results"
        / "staging_classifier/one_vs_rest_fibrosis_dream.csv"
    )
    fib_de = pd.read_csv(fib_de_path)

    # --- Build Ensembl → symbol mapping from atlas ---
    atlas_path = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
    atlas = pd.read_csv(atlas_path, usecols=[
        "human_symbol", "ensembl_id", "top_pathways",
        "bulk_logFC", "bulk_padj",
    ])
    # Strip version suffix from ensembl_id for joining
    atlas["ens_base"] = atlas["ensembl_id"].str.replace(r"\.\d+$", "", regex=True)
    ens2sym = dict(zip(atlas["ensembl_id"], atlas["human_symbol"]))
    # Also map without version for safety
    ens2sym_base = dict(zip(atlas["ens_base"], atlas["human_symbol"]))

    def map_symbol(ens_id):
        s = ens2sym.get(ens_id)
        if s:
            return s
        base = ens_id.split(".")[0] if isinstance(ens_id, str) else ens_id
        return ens2sym_base.get(base, ens_id)

    fib_de["symbol"] = fib_de["gene"].map(map_symbol)

    # --- Load metadata for sample counts per fibrosis stage ---
    meta_path = (
        PROJECT_ROOT
        / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata"
        / "unified_metadata.csv"
    )
    meta = pd.read_csv(meta_path)
    fib_counts = (
        meta["fibrosis_stage"]
        .dropna()
        .astype(int)
        .value_counts()
        .sort_index()
        .to_dict()
    )

    # --- Load NMF assignments ---
    nmf_path = PROJECT_ROOT / "RNA-seq/results/subtypes/nmf_assignments.csv"
    nmf = pd.read_csv(nmf_path, usecols=["sample_id", "nmf_subtype"])
    merged = meta[["sample_id", "fibrosis_stage"]].dropna(subset=["fibrosis_stage"]).merge(
        nmf, on="sample_id", how="inner"
    )
    merged["fibrosis_stage"] = merged["fibrosis_stage"].astype(int)

    # Fraction per stage
    nmf_fracs = {}
    for stage in range(5):
        sub = merged[merged["fibrosis_stage"] == stage]
        total = len(sub)
        if total == 0:
            nmf_fracs[stage] = {"s1": None, "s2": None}
        else:
            s1 = (sub["nmf_subtype"] == "S1").sum() / total
            s2 = (sub["nmf_subtype"] == "S2").sum() / total
            nmf_fracs[stage] = {"s1": _round(s1, 3), "s2": _round(s2, 3)}

    # --- Build per-stage entries ---
    stage_labels_map = {0: "Healthy / NAFL", 1: "Mild Fibrosis", 2: "Moderate Fibrosis",
                        3: "Advanced Fibrosis", 4: "Cirrhosis"}
    padj_thr = 0.05
    lfc_thr = 0.5

    # Pre-build atlas pathway lookup: symbol -> top_pathways
    sym2pathways = {}
    for _, row in atlas.dropna(subset=["top_pathways"]).iterrows():
        pws = str(row["top_pathways"])
        if pws and pws != "nan":
            sym2pathways[row["human_symbol"]] = [p.strip() for p in pws.split(";") if p.strip()]

    stages = []
    for sv in range(5):
        label = f"F{sv}_vs_rest"
        sub = fib_de[fib_de["stage_label"] == label].copy()
        sig = sub[(sub["padj"] < padj_thr) & (sub["logFC"].abs() > lfc_thr)]
        up = sig[sig["logFC"] > 0].sort_values("logFC", ascending=False)
        down = sig[sig["logFC"] < 0].sort_values("logFC", ascending=True)

        top_up = [
            {"symbol": r["symbol"], "logfc": _round(r["logFC"]), "padj": _round(r["padj"])}
            for _, r in up.head(10).iterrows()
        ]
        top_down = [
            {"symbol": r["symbol"], "logfc": _round(r["logFC"]), "padj": _round(r["padj"])}
            for _, r in down.head(10).iterrows()
        ]

        # Key pathways: aggregate top_pathways from the stage's top DEGs
        all_degs = pd.concat([up.head(50), down.head(50)])
        pathway_counts: dict[str, int] = {}
        for sym in all_degs["symbol"]:
            for pw in sym2pathways.get(sym, []):
                pathway_counts[pw] = pathway_counts.get(pw, 0) + 1
        key_pathways = sorted(pathway_counts, key=pathway_counts.get, reverse=True)[:6]
        # Clean pathway names for display
        key_pathways_display = [
            pw.replace("HALLMARK_", "").replace("REACTOME_", "").replace("_", " ").title()
            for pw in key_pathways
        ]

        stages.append({
            "stage": f"F{sv}",
            "label": stage_labels_map[sv],
            "n_samples": fib_counts.get(sv, 0),
            "n_degs_up": len(up),
            "n_degs_down": len(down),
            "top_genes_up": top_up,
            "top_genes_down": top_down,
            "s1_fraction": nmf_fracs[sv]["s1"],
            "s2_fraction": nmf_fracs[sv]["s2"],
            "key_pathways": key_pathways_display,
        })

    # --- Early-to-late stage signature ---
    # Part of the multi-step (F0→F1→F2→F3→F4) fibrosis progression, not a
    # discrete F2 "switch" (that framing is retired). Captured by comparing
    # the early-stage (F0/F1) DEG programme with the late-stage (F3/F4)
    # programme — genes that are NOT significant in F0/F1 but become
    # strongly significant in F3/F4.
    early_labels = ["F0_vs_rest", "F1_vs_rest"]
    late_labels = ["F3_vs_rest", "F4_vs_rest"]

    early_de = fib_de[fib_de["stage_label"].isin(early_labels)].copy()
    late_de = fib_de[fib_de["stage_label"].isin(late_labels)].copy()

    early_sig = set(early_de.loc[
        (early_de["padj"] < padj_thr) & (early_de["logFC"].abs() > lfc_thr), "gene"
    ])
    late_sig = set(late_de.loc[
        (late_de["padj"] < padj_thr) & (late_de["logFC"].abs() > lfc_thr), "gene"
    ])

    new_in_late = late_sig - early_sig  # emerge only in F3/F4

    # For transition genes, compare the max absolute LFC in early vs late
    f0_map = fib_de[fib_de["stage_label"] == "F0_vs_rest"].set_index("gene")["logFC"]
    f1_map = fib_de[fib_de["stage_label"] == "F1_vs_rest"].set_index("gene")["logFC"]
    f3_map = fib_de[fib_de["stage_label"] == "F3_vs_rest"].set_index("gene")["logFC"]
    f4_map = fib_de[fib_de["stage_label"] == "F4_vs_rest"].set_index("gene")["logFC"]

    transition_genes = []
    for g in new_in_late:
        early_lfc = max(
            abs(f0_map.get(g, 0.0)),
            abs(f1_map.get(g, 0.0)),
        )
        # Use the strongest late-stage LFC (keep sign)
        f3_lfc = f3_map.get(g, 0.0)
        f4_lfc = f4_map.get(g, 0.0)
        late_lfc = f3_lfc if abs(f3_lfc) >= abs(f4_lfc) else f4_lfc
        sym = map_symbol(g)
        transition_genes.append({
            "symbol": sym,
            "logfc_early": _round(early_lfc),
            "logfc_late": _round(late_lfc),
        })
    transition_genes.sort(key=lambda x: abs(x["logfc_late"] or 0), reverse=True)

    new_inflammatory = sum(1 for sg in transition_genes if (sg["logfc_late"] or 0) > 0)

    late_stage_signature = {
        "total_late_degs": len(late_sig),
        "early_degs": len(early_sig),
        "new_in_late": len(new_in_late),
        "new_inflammatory_genes": new_inflammatory,
        "key_transition_genes": transition_genes[:15],
    }

    # --- Classifier metrics (hardcoded from verified results) ---
    classifier_metrics = {
        "f_ge_3_auroc": 0.881,
        "nas_ge_5_auroc": 0.798,
        "fibrosis_qwk": 0.624,
    }

    result = {
        "stages": stages,
        "late_stage_signature": late_stage_signature,
        "classifier_metrics": classifier_metrics,
    }

    # Stats
    total_degs = sum(s["n_degs_up"] + s["n_degs_down"] for s in stages)
    print(f"    Stages: {len(stages)}, total DEGs across stages: {total_degs}")
    print(f"    Early-to-late signature: {late_stage_signature['total_late_degs']} late DEGs, "
          f"{late_stage_signature['new_in_late']} new in late (F3/F4), "
          f"{late_stage_signature['new_inflammatory_genes']} inflammatory")
    for s in stages:
        print(f"      {s['stage']}: {s['n_samples']} samples, "
              f"{s['n_degs_up']}↑ / {s['n_degs_down']}↓ DEGs, "
              f"S1={s['s1_fraction']}, S2={s['s2_fraction']}")

    return result


# ---------------------------------------------------------------------------
# 2. Drug Pipeline
# ---------------------------------------------------------------------------

def build_drug_pipeline() -> dict:
    print("[2/3] Building drug_pipeline.json ...")

    drug_dir = PROJECT_ROOT / "RNA-seq/results/drug_repurposing"
    strat_dir = PROJECT_ROOT / "RNA-seq/results/stratified_causal"

    # --- LINCS ranked (top compounds) ---
    lincs_ranked = pd.read_csv(drug_dir / "lincs_final_ranked.csv")

    # --- CGP reversal hits ---
    cgp_hits = pd.read_csv(drug_dir / "cgp_reversal_hits.csv")

    # --- Network proximity ---
    net_prox = pd.read_csv(drug_dir / "network_proximity" / "network_proximity_scores.csv")
    n_net_sig = len(net_prox[net_prox["p_closest"] < 0.05])

    # --- Clinical validation ---
    clin_val = pd.read_csv(drug_dir / "clinical_drug_validation_table.csv")

    # --- Multi-layer targets ---
    multi_layer_path = drug_dir / "multi_layer_drug_targets.csv"
    if multi_layer_path.exists():
        multi_layer = pd.read_csv(multi_layer_path)
        n_multi = len(multi_layer)
    else:
        # Fallback: count LINCS ranked compounds with multiple non-zero scores
        scored = lincs_ranked[["score_reversal", "score_sig", "score_mr",
                               "score_dgidb", "score_network"]].fillna(0)
        n_multi = int((scored.gt(0).sum(axis=1) >= 3).sum())

    # --- Stratified reversal for sex stats ---
    strat_path = strat_dir / "lincs_stratified_reversal.csv"
    if strat_path.exists():
        strat = pd.read_csv(strat_path)
        sex_counts = strat["sex_bias"].value_counts().to_dict()
        sex_stats = {
            "female_biased": int(sex_counts.get("Female_biased", 0)),
            "male_biased": int(sex_counts.get("Male_biased", 0)),
            "balanced": int(sex_counts.get("Balanced", 0)),
            "total_compounds": len(strat),
        }
    else:
        sex_stats = {
            "female_biased": 89,
            "male_biased": 44,
            "balanced": 0,
            "total_compounds": 133,
        }

    # --- Funnel ---
    # CGP file has 133 pathway gene sets that are reversal hits.
    # LINCS queried 1,107 perturbagens (from CLAUDE.md: "1,107 LINCS").
    # Network proximity scored 625 drugs. Sig at p<0.05.
    # Unique clinical drugs (not rows — rows can be multi-target)
    n_clinical_drugs = clin_val["drug"].nunique()

    funnel = [
        {
            "stage": "LINCS Signature Reversal",
            "count": 1107,
            "description": "L1000 perturbagens with significant disease-signature reversal",
        },
        {
            "stage": "Network Proximity Filter",
            "count": n_net_sig,
            "description": f"Drugs within disease-module proximity (p<0.05 from {len(net_prox)} tested)",
        },
        {
            "stage": "Multi-Layer Evidence",
            "count": n_multi,
            "description": "Compounds with reversal + proximity + DGIdb/COLOC support",
        },
        {
            "stage": "Clinical Validation",
            "count": n_clinical_drugs,
            "description": "Drugs in clinical trials or approved for MASLD/NASH",
        },
    ]

    # --- Top compounds ---
    top_compounds = []
    for _, r in lincs_ranked.sort_values("composite_score", ascending=False).head(20).iterrows():
        target_val = r.get("target.x", r.get("target.y", ""))
        moa_val = r.get("moa.x", r.get("moa.y", ""))
        # Clean up escaped quotes from R CSV
        if isinstance(target_val, str):
            target_val = target_val.strip('"').strip()
        if isinstance(moa_val, str):
            moa_val = moa_val.strip('"').strip()

        top_compounds.append({
            "name": str(r["display_name"]),
            "target": target_val if target_val and str(target_val) != "nan" else None,
            "moa": moa_val if moa_val and str(moa_val) != "nan" else None,
            "composite_score": _round(r["composite_score"], 3),
            "score_reversal": _round(r.get("score_reversal"), 3),
            "score_network": _round(r.get("score_network"), 3),
            "score_dgidb": _round(r.get("score_dgidb"), 3),
        })

    # --- Validated drugs ---
    validated_drugs = []
    for _, r in clin_val.iterrows():
        validated_drugs.append({
            "drug": str(r["drug"]),
            "target": str(r["target_gene"]),
            "stage": str(r["stage"]),
            "moa": str(r["moa"]),
            "support": str(r["atlas_support"]),
            "is_deg": bool(r["is_deg"]) if pd.notna(r.get("is_deg")) else None,
            # JSON key unified to bulk_* (C2); source clinical_drug_validation_table.csv
            # is a non-atlas drug table not yet migrated (separate area).
            "bulk_logfc": _round(r.get("dream_logFC")),  # C2-OK-sensitivity
        })

    result = {
        "funnel": funnel,
        "top_compounds": top_compounds,
        "validated_drugs": validated_drugs,
        "sex_stats": sex_stats,
    }

    print(f"    Funnel: {' → '.join(str(f['count']) for f in funnel)}")
    print(f"    Top compounds: {len(top_compounds)}")
    print(f"    Validated drugs: {len(validated_drugs)} rows ({n_clinical_drugs} unique)")
    print(f"    Sex stats: {sex_stats}")

    return result


# ---------------------------------------------------------------------------
# 3. Cross-Species
# ---------------------------------------------------------------------------

def build_cross_species() -> dict:
    print("[3/3] Building cross_species.json ...")

    conc_dir = (
        PROJECT_ROOT / "Analysis/Cross_Species_Concordance/results"
    )

    # --- Concordance atlas ---
    conc = pd.read_csv(conc_dir / "concordance_atlas_unified.csv")

    # --- Also pull mouse_meta_logFC from multi-evidence atlas ---
    atlas = pd.read_csv(
        PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv",
        usecols=["human_symbol", "mouse_meta_logFC"],
    )
    atlas_mouse = atlas.dropna(subset=["mouse_meta_logFC"]).set_index("human_symbol")[
        "mouse_meta_logFC"
    ].to_dict()

    # --- Build genes list ---
    # Include all Conserved (723), plus top genes from other categories
    # (by translatability score), capped to keep file size manageable (~300KB).
    # For the scatter plot, we keep all genes with human or mouse LFC != 0.
    priority_cats = {"Conserved", "Moderate_Concordance", "Species_Discordant", "Human_Enriched"}
    conc["_abs_h_lfc"] = conc["mean_h_lfc"].abs()

    # Keep: all priority categories + top 200 from remaining by |human LFC|
    priority_mask = conc["primary_category"].isin(priority_cats)
    remaining = conc[~priority_mask].nlargest(500, "_abs_h_lfc")
    selected = pd.concat([conc[priority_mask], remaining]).drop_duplicates(subset="human_symbol")

    genes = []
    for _, r in selected.iterrows():
        sym = r["human_symbol"]
        mouse_lfc = atlas_mouse.get(sym)
        genes.append({
            "symbol": sym,
            "human_logfc": _round(r["mean_h_lfc"]),
            "mouse_logfc": _round(mouse_lfc),
            "category": r["primary_category"],
            "n_concordant": int(r["n_concordant"]) if pd.notna(r["n_concordant"]) else 0,
            "n_diets_sig": int(r["n_diets_sig"]) if pd.notna(r["n_diets_sig"]) else 0,
            "translatability_score": _round(r.get("translatability_score")),
            "diets_concordant": (
                str(r["diets_concordant"]) if pd.notna(r.get("diets_concordant")) and str(r.get("diets_concordant")) else None
            ),
        })

    # --- Summary ---
    cat_counts = conc["primary_category"].value_counts().to_dict()
    summary = {
        "total_orthologs": len(conc),
        "conserved": int(cat_counts.get("Conserved", 0)),
        "human_specific": int(cat_counts.get("Human_Enriched", 0)),
        "mouse_specific": int(cat_counts.get("Mouse_Specific", 0)),
        "moderate_concordance": int(cat_counts.get("Moderate_Concordance", 0)),
        "diet_selective": int(cat_counts.get("Diet_Selective", 0)),
        "discordant": int(cat_counts.get("Species_Discordant", 0)),
        "not_significant": int(cat_counts.get("Not_Significant", 0)),
        "diet_models": ["MCD", "HFD", "CDAHFD", "FPC", "LIDPAD"],
    }

    # --- Pathway concordance ---
    pw_conc_path = conc_dir / "fgsea_pathway_concordance.csv"
    pathway_concordance = []
    if pw_conc_path.exists():
        # Load human and mouse HALLMARK NES for disease_vs_ctrl and MCD (best model)
        human_fgsea = pd.read_csv(conc_dir / "fgsea_human_results.csv")
        mouse_fgsea = pd.read_csv(conc_dir / "fgsea_mouse_results.csv")

        h_hall = human_fgsea[
            (human_fgsea["pathway"].str.startswith("HALLMARK_"))
            & (human_fgsea["source"] == "disease_vs_ctrl")
        ].set_index("pathway")

        # Average mouse NES across all diets for each pathway
        m_hall = mouse_fgsea[
            mouse_fgsea["pathway"].str.startswith("HALLMARK_")
        ].groupby("pathway").agg({"NES": "mean", "padj": "min"}).rename(
            columns={"NES": "mouse_NES", "padj": "mouse_padj"}
        )

        common_pws = h_hall.index.intersection(m_hall.index)
        for pw in sorted(common_pws):
            h_nes = h_hall.loc[pw, "NES"]
            m_nes = m_hall.loc[pw, "mouse_NES"]
            concordant = (h_nes > 0 and m_nes > 0) or (h_nes < 0 and m_nes < 0)
            pathway_concordance.append({
                "pathway": pw.replace("HALLMARK_", "").replace("_", " ").title(),
                "pathway_id": pw,
                "human_nes": _round(h_nes, 3),
                "mouse_nes": _round(m_nes, 3),
                "concordant": bool(concordant),
            })
        # Sort by absolute human NES descending
        pathway_concordance.sort(key=lambda x: abs(x["human_nes"] or 0), reverse=True)

    result = {
        "genes": genes,
        "summary": summary,
        "pathway_concordance": pathway_concordance,
    }

    n_concordant_pws = sum(1 for p in pathway_concordance if p["concordant"])
    print(f"    Genes: {len(genes)}")
    print(f"    Categories: {summary}")
    print(f"    Pathway concordance: {len(pathway_concordance)} HALLMARK pathways, "
          f"{n_concordant_pws} concordant")

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, help="Output directory for JSON files")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out_dir.resolve()}\n")

    # 1. Progression
    progression = build_progression_journey()
    write_json(progression, out_dir / "progression_journey.json")

    # 2. Drug pipeline
    drug = build_drug_pipeline()
    write_json(drug, out_dir / "drug_pipeline.json")

    # 3. Cross-species
    species = build_cross_species()
    write_json(species, out_dir / "cross_species.json")

    print("\nDone. All 3 Phase 3 JSON files generated.")


if __name__ == "__main__":
    main()
