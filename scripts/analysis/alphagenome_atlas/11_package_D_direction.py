#!/usr/bin/env python3
"""Step 11 (Package D): directional consistency and disagreement, three separate comparisons.

1. Predicted allele effect (Atlas liver RNA, signed, ref→alt) vs measured eQTL effect (Broadaway,
   harmonised by the continuum-translation rules: allele_state, aligned_eqtl_beta). Per variant on the
   signal's target gene; posterior-weighted concordance per signal; locus-balanced by analysis block.
   Also the prespecified replication of the archived 211-variant AlphaGenome gate (65b), same set,
   same orientation rule, same statistics, seed 42. Written prediction: at chance.
2. Risk-allele-associated expression direction (marginal_expression_direction from the translation
   build, resolved only at ≥0.95 mass) vs disease-associated expression (continuum translation
   remodeling_state / axis betas). Reported per signal; unresolved kept.
3. Bulk disease effect vs within-lineage disease effect: deferred to the existing atac-context
   program tables (not re-derived here); recorded as not_applicable in this step.

Outputs (tables/): direction_variant_level.tsv.gz, direction_comparisons.tsv, ag_atlas_gate_replication.json
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict

import numpy as np
import pandas as pd

import lib_atlas as la

P = la.prespec()
TABLES = la.out_root() / "tables"
RAW = la.out_root() / "raw"
TR = la.PROJECT / "RNA-seq/results/histology_anchored_continuum/translation/hac-translation-20260906T230500Z"
AG_SCORES = la.PROJECT / "GWAS/finemapping/results/seqfunc/alphagenome_eqtl_scores.tsv"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
TRUTH = la.PROJECT / "GWAS/finemapping/results/seqfunc/broadaway_benchmark_truth.tsv"
CUT = P["coverage"]["direction_cutoff"]
GATE_SEED = P["seeds"]["gate_replication"]


def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    order = np.argsort(scores)
    ranks = np.empty(len(scores)); ranks[order] = np.arange(1, len(scores) + 1)
    # average ranks for ties
    s_sorted = scores[order]
    i = 0
    while i < len(s_sorted):
        j = i
        while j + 1 < len(s_sorted) and s_sorted[j + 1] == s_sorted[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2
        i = j + 1
    pos = labels == 1
    n1, n0 = pos.sum(), (~pos).sum()
    return float((ranks[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main() -> None:
    signals = {s["signal_uid"]: s for s in la.read_tsv(TABLES / "eligible_signals.tsv")}
    liver = pd.read_csv(TABLES / "liver_summaries.tsv.gz", sep="\t", keep_default_na=False)
    rna = liver[liver["scorer"] == "RNA_SEQ"].set_index(["variant_uid", "gene_id"])
    xw = pd.read_csv(TABLES / "variant_crosswalk.tsv.gz", sep="\t", keep_default_na=False).set_index("source_variant_id")

    # --- comparison 1: per-variant predicted vs measured eQTL direction, signals of universes A/C
    var_rows = []
    per_signal = defaultdict(lambda: {"pos": 0.0, "neg": 0.0, "unres": 0.0, "n": 0, "n_scored": 0})
    study_cache: dict[str, pd.DataFrame] = {}
    for sig, s in signals.items():
        if s["universe"] == "B_direct":
            continue
        path = TR / "signed_variants" / f"{s['gwas_name']}.tsv.gz"
        if not path.exists():
            continue
        if s["gwas_name"] not in study_cache:
            df = pd.read_csv(path, sep="\t", keep_default_na=False)
            df["_key"] = df["ensembl"].astype(str) + "|" + df["signal_pair_index"].astype(str)
            study_cache[s["gwas_name"]] = {k: g for k, g in df.groupby("_key")}
        sv = study_cache[s["gwas_name"]].get(f"{s['ensembl']}|{s['signal_pair_index']}", pd.DataFrame())
        for _, r in sv.iterrows():
            src = f"hg19:{r['chr']}:{r['position']}:{r['allele1']}:{r['allele2']}"
            w = float(r["SNP.PP.H4"])
            entry = per_signal[sig]; entry["n"] += 1
            if src not in xw.index or xw.at[src, "mapping_status"] != "mapped":
                entry["unres"] += w
                continue
            uid = xw.at[src, "variant_uid"]
            swap = xw.at[src, "allele_swap"]
            key = (uid, s["ensembl"])
            if key not in rna.index or r["allele_state"] not in ("direct", "swapped"):
                entry["unres"] += w
                var_rows.append({"signal_uid": sig, "variant_uid": uid, "weight": w, "allele_state": r["allele_state"], "atlas_scored": key in rna.index,
                                 "predicted_direction": "", "measured_direction": "", "concordant": ""})
                continue
            # Atlas: ref→alt on hg38. Translation: aligned_eqtl_beta is oriented to allele1 (GWAS effect allele) after harmonisation.
            # crosswalk allele_swap says whether allele1 is hg38_alt ('none' → allele1==ref? no: allele_swap 'none' means allele1 == ref).
            pred = float(rna.at[key, "liver_median_raw"]) if not isinstance(rna.at[key, "liver_median_raw"], pd.Series) else float(rna.at[key, "liver_median_raw"].iloc[0])
            pred_a1 = -pred if swap == "none" else pred   # effect of allele1 relative to allele2
            eqtl_a1 = float(r["aligned_eqtl_beta"])
            pd_, md = np.sign(pred_a1), np.sign(eqtl_a1)
            conc = bool(pd_ == md) if pd_ != 0 and md != 0 else None
            if conc is None:
                entry["unres"] += w
            elif conc:
                entry["pos"] += w
            else:
                entry["neg"] += w
            entry["n_scored"] += 1
            var_rows.append({"signal_uid": sig, "variant_uid": uid, "weight": w, "allele_state": r["allele_state"], "atlas_scored": True,
                             "predicted_direction": int(pd_), "measured_direction": int(md), "concordant": conc,
                             "pred_allele1_liver_rna": pred_a1, "eqtl_beta_allele1": eqtl_a1, "gwas_beta_allele1": r["gwas_beta"]})
    la.write_tsv_once(TABLES / "direction_variant_level.tsv.gz", var_rows, ["signal_uid", "variant_uid", "weight", "allele_state", "atlas_scored",
                      "predicted_direction", "measured_direction", "concordant", "pred_allele1_liver_rna", "eqtl_beta_allele1", "gwas_beta_allele1"])

    # --- comparison 2 inputs: risk-allele expression direction (translation) and disease remodeling state (translation catalog)
    cat = pd.read_csv(TR / "gene_catalog_translation.tsv.gz", sep="\t", keep_default_na=False).drop_duplicates("ensembl").set_index("ensembl")
    comp_rows = []
    for sig, s in signals.items():
        if s["universe"] == "B_direct":
            comp_rows.append({"signal_uid": sig, "universe": s["universe"], "analysis_block": s["analysis_block"], "comparison": "predicted_vs_measured_eqtl",
                              "state": "not_applicable_no_colocalized_gene"})
            continue
        e = per_signal[sig]
        tot = e["pos"] + e["neg"] + e["unres"]
        resolved = tot > 0 and max(e["pos"], e["neg"]) >= CUT * tot
        comp_rows.append({"signal_uid": sig, "universe": s["universe"], "analysis_block": s["analysis_block"], "comparison": "predicted_vs_measured_eqtl",
                          "n_variants": e["n"], "n_scored": e["n_scored"], "concordant_mass": e["pos"], "discordant_mass": e["neg"], "unresolved_mass": e["unres"],
                          "state": ("concordant" if e["pos"] >= e["neg"] else "discordant") if resolved else "unresolved",
                          "rule": f"posterior mass ≥ {CUT} of original total in one direction; no renormalisation; palindromes unresolved"})
        # comparison 2
        g = cat.loc[s["ensembl"]] if s["ensembl"] in cat.index else None
        row = {"signal_uid": sig, "universe": s["universe"], "analysis_block": s["analysis_block"], "comparison": "risk_expression_vs_disease_expression"}
        if g is None:
            row["state"] = "untested_gene_absent_from_translation_catalog"
        else:
            row.update({"genetic_evidence_state": g["genetic_evidence_state"], "remodeling_state": g["remodeling_state"],
                        "continuum_associated": g["continuum_associated"], "beta_fixed_projection": g["beta_fixed_projection"],
                        "protein_state": g["protein_state"], "target_readout_identity": g["target_readout_identity"],
                        "state": "reported_from_translation_build_no_new_test", "source": str(TR.relative_to(la.PROJECT))})
        comp_rows.append(row)
        comp_rows.append({"signal_uid": sig, "universe": s["universe"], "analysis_block": s["analysis_block"], "comparison": "bulk_vs_within_lineage_disease_effect",
                          "state": "not_rederived_here; see atac-context-v3 program_ATAC_native and program_meta_effects (translation build)"})
    la.write_tsv_once(TABLES / "direction_comparisons.tsv", comp_rows, sorted({k for r in comp_rows for k in r}))

    # --- gate replication on the archived 211-variant set with Atlas liver RNA scores
    import pysam
    fasta = pysam.FastaFile(FASTA)
    truth = la.gate_set(la.read_tsv(TRUTH), la.read_tsv(AG_SCORES), fasta.fetch)
    gate_dir = RAW / "atlas_gate"
    gate = {}
    if gate_dir.exists():
        import anndata
        import atlas_archive as aa
        for c in sorted(gate_dir.glob("chunk_*")):
            a = anndata.read_h5ad(c / "RNA_SEQ.h5ad")
            liver_mask = a.var["ontology_curie"].isin(["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]).values
            obs = a.obs.reset_index(drop=True)
            x = np.asarray(a.X)[:, liver_mask]
            for i in range(len(obs)):
                gate[(aa.variant_uid_from_str(str(obs.at[i, "variant"])), str(obs.at[i, "gene_id"]).split(".")[0])] = float(np.nanmean(x[i]))
    scores, labels, n_unscored = [], [], 0
    for r in truth:
        v = gate.get((r["variant_uid"], r["ensembl"]))
        if v is None or math.isnan(v):
            n_unscored += 1
            continue
        scores.append(r["orientation_sign"] * v); labels.append(r["label"])
    result = {"set": "broadaway_benchmark_truth is_signal_lead & !strand_ambiguous", "n_leads": len(truth), "n_scored": len(scores), "n_unscored": n_unscored,
              "prior_verdict": {"auroc": 0.5607, "p": 0.066, "tier": "FAIL", "n": 211}, "written_prediction": "at chance",
              "note": "Atlas liver RNA_SEQ mean over primary-liver+hepatocyte tracks, hg38 ref→alt from the 65-series score file (FASTA-validated), oriented to the eQTL effect allele exactly as 65b (effect==hg38_alt → +, ==hg38_ref → −); SNV leads only"}
    if len(scores) >= 20:
        s_, l_ = np.array(scores), np.array(labels)
        rng = np.random.default_rng(GATE_SEED)
        obs = auroc(s_, l_)
        boots = [auroc(s_[i], l_[i]) for i in (rng.integers(0, len(s_), len(s_)) for _ in range(2000))]
        perms = np.array([auroc(s_, rng.permutation(l_)) for _ in range(10000)])
        result.update({"auroc": obs, "ci95": [float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))],
                       "perm_p": float(((perms >= obs).sum() + 1) / (len(perms) + 1)), "sign_concordance": float(((s_ > 0) == (l_ == 1)).mean()),
                       "label_balance": {"pos": int(l_.sum()), "neg": int((l_ == 0).sum())}, "seed": GATE_SEED})
    else:
        result["state"] = "gate archive absent or too small; run 06_query_atlas.py gate"
    json.dump(result, (TABLES / "ag_atlas_gate_replication.json").open("w"), indent=1)
    la.log(f"Package D: {len(var_rows)} variant rows, {len(comp_rows)} comparison rows; gate n_scored={len(scores)} auroc={result.get('auroc')}")


if __name__ == "__main__":
    main()
