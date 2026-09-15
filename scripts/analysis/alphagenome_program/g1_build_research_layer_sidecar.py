#!/usr/bin/env python
"""G. Integration: build the candidate research-layer sidecar for the Gene Catalog.

Assembly only. Nothing here recomputes a statistic; every value is read from a deposited
table under GWAS/finemapping/results/alphagenome_program/ and the source file is recorded
in the column dictionary. Every row is marked candidate_not_adopted.

Key rule (spec section 8, wave-5 G assignment): key on (ensembl_id, credible_set_id) for
signals that carry a gene, and on signal_uid otherwise. The 173 B_direct SuSiE signals
carry no ensembl gene, so they take the signal_uid key.

Measurements and model predictions are kept in separate columns. Every signed column
names the allele or contrast its sign refers to.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PROG = ROOT / "GWAS/finemapping/results/alphagenome_program"
B1 = PROG / "b1-response-layer-v2-20260914T204700Z"
OUT = PROG / "g-integration-20260915T072500Z"

SIGNALS = B1 / "response_layer_signals.tsv"
GATE = B1 / "coverage_gate.tsv"
VARIANTS = B1 / "response_layer_variants.tsv.gz"

MEAS_PRED_NOTE = (
    "meas_* columns are measured quantities; pred_* columns are model predictions. "
    "They are never pooled and are not on one scale."
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def max_abs(series: pd.Series) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return np.nan
    return float(s.loc[s.abs().idxmax()])


def signed_max_abs(frame: pd.DataFrame, col: str) -> pd.Series:
    """Per signal_uid, the value with the largest absolute magnitude (sign preserved)."""
    sub = frame[["signal_uid", col]].copy()
    sub[col] = pd.to_numeric(sub[col], errors="coerce")
    sub = sub.dropna(subset=[col])
    if sub.empty:
        return pd.Series(dtype=float)
    return sub.groupby("signal_uid")[col].apply(lambda s: float(s.loc[s.abs().idxmax()]))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    sig = pd.read_csv(SIGNALS, sep="\t", dtype=str)
    gate = pd.read_csv(GATE, sep="\t", dtype=str)
    var = pd.read_csv(VARIANTS, sep="\t", dtype=str, low_memory=False)

    assert len(sig) == 989, f"expected 989 signals, got {len(sig)}"
    assert len(gate) == 989
    assert len(var) == 10511, f"expected 10511 variant rows, got {len(var)}"
    assert set(sig.atlas_extension_state.unique()) == {"candidate_not_adopted"}

    # ---- sign-convention strings, read from the layer, never invented -------------
    def one(col: str) -> str:
        vals = var[col].dropna().unique()
        assert len(vals) == 1, f"{col} carries {len(vals)} conventions"
        return str(vals[0])

    ase_sign = one("ase_sign_refers_to")
    api_sign = one("model_api_sign_refers_to")
    dir_sign = one("direction_sign_refers_to")
    da_sign_67 = one("measured_atac_logFC_gse281367_sign_refers_to")
    da_sign_32 = one("measured_atac_logFC_gse244832_sign_refers_to")

    # ---- per-signal aggregates from the variant table ----------------------------
    agg = pd.DataFrame(index=pd.Index(sig.signal_uid, name="signal_uid"))

    meas_cols = {
        "meas_allelic_imbalance_gse281367_max_abs_log2_alt_over_ref": "ase_gse281367_mean_log2_alt_over_ref",
        "meas_allelic_imbalance_gse244832_max_abs_log2_alt_over_ref": "ase_gse244832_mean_log2_alt_over_ref",
        "meas_atac_da_gse281367_max_abs_log2FC": "measured_atac_logFC_gse281367",
        "meas_atac_da_gse244832_max_abs_log2FC": "measured_atac_logFC_gse244832",
        "meas_liver_eqtl_max_abs_beta_allele1": "eqtl_beta_allele1",
        "meas_gwas_max_abs_beta_allele1": "gwas_beta_allele1",
    }
    pred_cols = {
        "pred_atlas_atac_primary_liver_max_abs_quantile": "atac_primary_liver_quantile",
        "pred_atlas_dnase_primary_liver_max_abs_quantile": "dnase_primary_liver_quantile",
        "pred_atlas_h3k27ac_primary_liver_max_abs_quantile": "h3k27ac_primary_liver_quantile",
        "pred_atlas_rna_target_gene_primary_liver_max_abs_quantile": "rna_target_gene_primary_liver_quantile",
        "pred_atlas_splice_site_usage_primary_liver_max_abs_quantile": "splice_site_usage_primary_liver_quantile",
        "pred_model_api_rna_max_abs_log2": "model_api_rna_log2",
        "pred_model_api_splice_max_abs_log2": "model_api_splice_log2",
        "pred_model_api_atac_max_abs_log2": "model_api_atac_log2",
        "pred_model_api_dnase_max_abs_log2": "model_api_dnase_log2",
        "pred_model_api_h3k27ac_max_abs_log2": "model_api_h3k27ac_log2",
    }
    for out_col, src in {**meas_cols, **pred_cols}.items():
        agg[out_col] = signed_max_abs(var, src)

    # het-donor counts behind the allelic-imbalance measurement (donors are the unit)
    for out_col, src in [
        ("meas_allelic_imbalance_gse281367_max_n_het_donors", "ase_gse281367_n_het_donors"),
        ("meas_allelic_imbalance_gse244832_max_n_het_donors", "ase_gse244832_n_het_donors"),
    ]:
        s = pd.to_numeric(var[src], errors="coerce")
        agg[out_col] = pd.DataFrame({"signal_uid": var.signal_uid, "v": s}).dropna().groupby("signal_uid")["v"].max()

    da_state = var[["signal_uid", "measured_da_evidence_state"]].fillna("")
    da_state = da_state[da_state.measured_da_evidence_state.str.strip() != ""]
    agg["meas_da_evidence_states"] = da_state.groupby("signal_uid")["measured_da_evidence_state"].apply(
        lambda s: ";".join(sorted(set(";".join(s).split(";"))))
    )

    lineages = var[["signal_uid", "measured_lineages"]].fillna("")
    lineages = lineages[lineages.measured_lineages.str.strip() != ""]
    agg["meas_donor_snatac_peak_lineages"] = lineages.groupby("signal_uid")["measured_lineages"].apply(
        lambda s: ";".join(sorted(set(";".join(s).split(";"))))
    )

    agg = agg.reindex(sig.signal_uid.values)

    # ---- assemble ---------------------------------------------------------------
    g = gate.set_index("signal_uid")
    s = sig.set_index("signal_uid")

    out = pd.DataFrame(index=s.index)
    out["ensembl_id"] = s["ensembl"].fillna("")
    out["credible_set_id"] = out.index
    out["signal_uid"] = out.index
    out["sidecar_key"] = np.where(
        out["ensembl_id"].str.strip() != "",
        out["ensembl_id"].astype(str) + "|" + out["credible_set_id"].astype(str),
        out["signal_uid"].astype(str),
    )
    out["sidecar_key_type"] = np.where(
        out["ensembl_id"].str.strip() != "", "ensembl_id_and_credible_set_id", "signal_uid_only"
    )
    out["gene_symbol"] = s["gene"].fillna("")
    out["universe"] = s["universe"]
    out["trait_class"] = s["trait_class"]
    out["trait"] = s["trait"]
    out["gwas_name"] = s["gwas_name"]
    out["posterior_definition"] = s["posterior_definition"]
    out["analysis_block"] = s["analysis_block"].fillna("")
    out["n_distinct_1mb_blocks"] = s["n_distinct_1mb_blocks"]

    # release state, on every row
    out["release_state"] = "candidate_not_adopted"
    out["atlas_extension_state"] = s["atlas_extension_state"]
    out["changes_any_adopted_resource_number"] = "no"

    # coverage state
    out["atlas_coverage_state"] = s["atlas_coverage_state"]
    out["signal_gate"] = g["signal_gate"]
    out["queried_share"] = s["queried_share"]
    out["total_mass"] = s["total_mass"]
    out["excluded_share_indel_deferred"] = g["excluded_share_indel_deferred"]
    out["excluded_share_liftover_failed"] = g["excluded_share_liftover_failed"]
    out["top_variant_served"] = g["top_variant_served"]
    out["n_exported_variant_rows"] = s["n_exported_rows"]
    out["n_variant_rows_atlas_served"] = s["n_rows_served"]
    out["n_variant_rows_readable"] = s["n_rows_readable"]

    # distinguishable explanation classes
    out["classes_distinguishable_from_measurement"] = s["classes_distinguishable_from_measurement"]
    out["n_classes_distinguishable_from_measurement"] = s["n_classes_distinguishable_from_measurement"]
    out["class_coding_state"] = s["class_coding"]
    out["class_local_regulatory_state"] = s["class_local_regulatory"]
    out["class_rna_processing_state"] = s["class_rna_processing"]
    out["class_candidate_long_range_state"] = s["class_candidate_long_range"]
    out["locus_state"] = s["locus_state"]

    # measured channels: counts (from the B1 signals table) and values (aggregated here)
    out["meas_n_variants_in_donor_snatac_peak"] = s["n_rows_in_measured_peak"]
    out["meas_n_variants_with_da_evidence"] = s["n_rows_with_da_evidence"]
    out["meas_n_variants_with_allelic_imbalance"] = s["n_rows_with_measured_allelic_imbalance"]
    out["meas_n_variants_with_liver_eqtl_direction"] = s["n_rows_with_eqtl_direction"]
    for c in [
        "meas_allelic_imbalance_gse281367_max_abs_log2_alt_over_ref",
        "meas_allelic_imbalance_gse281367_max_n_het_donors",
        "meas_allelic_imbalance_gse244832_max_abs_log2_alt_over_ref",
        "meas_allelic_imbalance_gse244832_max_n_het_donors",
        "meas_atac_da_gse281367_max_abs_log2FC",
        "meas_atac_da_gse244832_max_abs_log2FC",
        "meas_da_evidence_states",
        "meas_donor_snatac_peak_lineages",
        "meas_liver_eqtl_max_abs_beta_allele1",
        "meas_gwas_max_abs_beta_allele1",
    ]:
        out[c] = agg[c].values

    counts = pd.DataFrame(
        {
            "peak": pd.to_numeric(s["n_rows_in_measured_peak"], errors="coerce").fillna(0),
            "da": pd.to_numeric(s["n_rows_with_da_evidence"], errors="coerce").fillna(0),
            "ase": pd.to_numeric(s["n_rows_with_measured_allelic_imbalance"], errors="coerce").fillna(0),
            "eqtl": pd.to_numeric(s["n_rows_with_eqtl_direction"], errors="coerce").fillna(0),
        }
    )
    names = {
        "peak": "donor_snatac_peak_overlap",
        "da": "donor_differential_accessibility_state",
        "ase": "measured_allelic_imbalance",
        "eqtl": "liver_eqtl_direction",
    }
    out["meas_any_measured_channel"] = np.where(counts.sum(axis=1) > 0, "yes", "no")
    out["meas_channels_present"] = [
        ";".join(names[k] for k in counts.columns if row[k] > 0) or "none"
        for _, row in counts.iterrows()
    ]

    # sign conventions of the measured columns
    out["meas_allelic_imbalance_sign_refers_to"] = np.where(
        out["meas_allelic_imbalance_gse281367_max_abs_log2_alt_over_ref"].notna()
        | out["meas_allelic_imbalance_gse244832_max_abs_log2_alt_over_ref"].notna(),
        ase_sign,
        "",
    )
    out["meas_atac_da_gse281367_sign_refers_to"] = np.where(
        out["meas_atac_da_gse281367_max_abs_log2FC"].notna(), da_sign_67, ""
    )
    out["meas_atac_da_gse244832_sign_refers_to"] = np.where(
        out["meas_atac_da_gse244832_max_abs_log2FC"].notna(), da_sign_32, ""
    )
    out["meas_eqtl_and_gwas_sign_refers_to"] = np.where(
        out["meas_liver_eqtl_max_abs_beta_allele1"].notna() | out["meas_gwas_max_abs_beta_allele1"].notna(),
        dir_sign,
        "",
    )

    # model predictions
    for c in pred_cols:
        out[c] = agg[c].values
    out["pred_n_variants_with_model_api_rescue"] = s["n_rows_with_model_api_rescue"]
    out["pred_atlas_quantile_sign_refers_to"] = np.where(
        out[[c for c in pred_cols if c.startswith("pred_atlas") and "quantile" in c]].notna().any(axis=1),
        "signed Atlas quantile of the predicted allele effect, ALT relative to REF in hg38 "
        "reference orientation; a quantile ranks a change against a distribution of changes "
        "and is not an effect size",
        "",
    )
    out["pred_model_api_sign_refers_to"] = np.where(
        out[[c for c in pred_cols if c.startswith("pred_model_api")]].notna().any(axis=1), api_sign, ""
    )
    out["pred_atlas_liver_panel_note"] = (
        "ATAC 3 adult tracks; DNase 2 adult + 1 embryonic; H3K27ac 2 adult; RNA and splice-site "
        "usage 3 adult + 1 child + 1 embryonic. Panel composition from "
        "b1-response-layer-v2-20260914T204700Z/track_panel_composition.tsv"
    )
    out["pred_next_experiment_rule_id"] = s["next_experiment_rule_id"]
    out["pred_candidate_mechanism_atlas_predicted"] = s["candidate_mechanism_atlas_predicted"]

    out["source_deposit"] = "b1-response-layer-v2-20260914T204700Z"
    out["measurement_vs_prediction_note"] = MEAS_PRED_NOTE

    out = out.reset_index(drop=True)

    # ---- guards ------------------------------------------------------------------
    assert len(out) == 989
    assert (out.release_state == "candidate_not_adopted").all()
    assert out.sidecar_key.is_unique, "sidecar key is not unique"
    n_gene = int((out.sidecar_key_type == "ensembl_id_and_credible_set_id").sum())
    n_uid = int((out.sidecar_key_type == "signal_uid_only").sum())
    assert n_gene + n_uid == 989
    assert n_uid == int((out.ensembl_id.str.strip() == "").sum())

    sidecar_path = OUT / "research_layer_sidecar.tsv"
    out.to_csv(sidecar_path, sep="\t", index=False)

    # ---- data dictionary ---------------------------------------------------------
    dict_rows = []
    kind_of = {}
    for c in out.columns:
        if c.startswith("meas_"):
            kind_of[c] = "measurement"
        elif c.startswith("pred_"):
            kind_of[c] = "model_prediction"
        else:
            kind_of[c] = "identity_or_state"
    units = {
        "meas_allelic_imbalance_gse281367_max_abs_log2_alt_over_ref": "log2(ALT/REF), mean over het donors",
        "meas_allelic_imbalance_gse244832_max_abs_log2_alt_over_ref": "log2(ALT/REF), mean over het donors",
        "meas_atac_da_gse281367_max_abs_log2FC": "log2 fold change, disease state contrast",
        "meas_atac_da_gse244832_max_abs_log2FC": "log2 fold change, disease state contrast",
        "meas_liver_eqtl_max_abs_beta_allele1": "eQTL beta",
        "meas_gwas_max_abs_beta_allele1": "trait-specific GWAS beta",
        "pred_model_api_rna_max_abs_log2": "log2(ALT/REF), hosted model API",
        "pred_model_api_splice_max_abs_log2": "log2(ALT/REF), hosted model API",
        "pred_model_api_atac_max_abs_log2": "log2(ALT/REF), hosted model API",
        "pred_model_api_dnase_max_abs_log2": "log2(ALT/REF), hosted model API",
        "pred_model_api_h3k27ac_max_abs_log2": "log2(ALT/REF), hosted model API",
        "pred_atlas_atac_primary_liver_max_abs_quantile": "signed Atlas quantile, dimensionless",
        "pred_atlas_dnase_primary_liver_max_abs_quantile": "signed Atlas quantile, dimensionless",
        "pred_atlas_h3k27ac_primary_liver_max_abs_quantile": "signed Atlas quantile, dimensionless",
        "pred_atlas_rna_target_gene_primary_liver_max_abs_quantile": "signed Atlas quantile, dimensionless",
        "pred_atlas_splice_site_usage_primary_liver_max_abs_quantile": "signed Atlas quantile, dimensionless",
    }
    sign_of = {
        "meas_allelic_imbalance_gse281367_max_abs_log2_alt_over_ref": "meas_allelic_imbalance_sign_refers_to",
        "meas_allelic_imbalance_gse244832_max_abs_log2_alt_over_ref": "meas_allelic_imbalance_sign_refers_to",
        "meas_atac_da_gse281367_max_abs_log2FC": "meas_atac_da_gse281367_sign_refers_to",
        "meas_atac_da_gse244832_max_abs_log2FC": "meas_atac_da_gse244832_sign_refers_to",
        "meas_liver_eqtl_max_abs_beta_allele1": "meas_eqtl_and_gwas_sign_refers_to",
        "meas_gwas_max_abs_beta_allele1": "meas_eqtl_and_gwas_sign_refers_to",
        "pred_model_api_rna_max_abs_log2": "pred_model_api_sign_refers_to",
        "pred_model_api_splice_max_abs_log2": "pred_model_api_sign_refers_to",
        "pred_model_api_atac_max_abs_log2": "pred_model_api_sign_refers_to",
        "pred_model_api_dnase_max_abs_log2": "pred_model_api_sign_refers_to",
        "pred_model_api_h3k27ac_max_abs_log2": "pred_model_api_sign_refers_to",
        "pred_atlas_atac_primary_liver_max_abs_quantile": "pred_atlas_quantile_sign_refers_to",
        "pred_atlas_dnase_primary_liver_max_abs_quantile": "pred_atlas_quantile_sign_refers_to",
        "pred_atlas_h3k27ac_primary_liver_max_abs_quantile": "pred_atlas_quantile_sign_refers_to",
        "pred_atlas_rna_target_gene_primary_liver_max_abs_quantile": "pred_atlas_quantile_sign_refers_to",
        "pred_atlas_splice_site_usage_primary_liver_max_abs_quantile": "pred_atlas_quantile_sign_refers_to",
    }
    src_of = {c: "b1-response-layer-v2-20260914T204700Z/response_layer_variants.tsv.gz" for c in units}
    for c in out.columns:
        dict_rows.append(
            {
                "column": c,
                "kind": kind_of[c],
                "unit": units.get(c, ""),
                "sign_convention_column": sign_of.get(c, ""),
                "source_file": src_of.get(c, "b1-response-layer-v2-20260914T204700Z/response_layer_signals.tsv"
                                          if c in set(s.columns) | {"signal_uid"} else
                                          "b1-response-layer-v2-20260914T204700Z/coverage_gate.tsv"
                                          if c in set(g.columns) else "assembled in this package"),
                "n_non_null": int(out[c].notna().sum() - (out[c].astype(str).str.strip() == "").sum()),
            }
        )
    dict_path = OUT / "research_layer_sidecar_dictionary.tsv"
    pd.DataFrame(dict_rows).to_csv(dict_path, sep="\t", index=False)

    summary = {
        "rows": int(len(out)),
        "key_ensembl_and_credible_set": n_gene,
        "key_signal_uid_only": n_uid,
        "universe": out.universe.value_counts().to_dict(),
        "release_state_unique": sorted(out.release_state.unique().tolist()),
        "atlas_coverage_state": out.atlas_coverage_state.value_counts().to_dict(),
        "any_measured_channel": out.meas_any_measured_channel.value_counts().to_dict(),
        "channels_present": out.meas_channels_present.value_counts().to_dict(),
        "n_classes_distinguishable": out.n_classes_distinguishable_from_measurement.value_counts().to_dict(),
        "signals_with_measured_allelic_imbalance_value": int(
            (
                out.meas_allelic_imbalance_gse281367_max_abs_log2_alt_over_ref.notna()
                | out.meas_allelic_imbalance_gse244832_max_abs_log2_alt_over_ref.notna()
            ).sum()
        ),
        "signals_with_model_api_rescue": int(
            pd.to_numeric(out.pred_n_variants_with_model_api_rescue, errors="coerce").fillna(0).gt(0).sum()
        ),
        "sign_conventions": {
            "measured_allelic_imbalance": ase_sign,
            "measured_atac_da_gse281367": da_sign_67,
            "measured_atac_da_gse244832": da_sign_32,
            "measured_eqtl_and_gwas": dir_sign,
            "model_api": api_sign,
        },
        "sidecar_sha256": sha256(sidecar_path),
        "dictionary_sha256": sha256(dict_path),
    }
    (OUT / "research_layer_sidecar_summary.json").write_text(json.dumps(summary, indent=2))

    # ---- MANIFEST of inputs ------------------------------------------------------
    inputs = [
        SIGNALS,
        GATE,
        VARIANTS,
        B1 / "response_layer_summary.json",
        B1 / "column_dictionary.tsv",
        B1 / "track_panel_composition.tsv",
        PROG / "a1-mpra-reliability-v2-20260914T204500Z/tables/reliability_by_context.tsv",
        PROG / "a1-mpra-reliability-v2-20260914T204500Z/tables/decision_rule_section3.tsv",
        PROG / "a2-rights-exposure-20260914T193626Z/tables/rights_exposure.tsv",
        PROG / "b3-saturation-20260914T193641Z/tables/b3_contrasts.tsv",
        PROG / "c1-matched-endogenous-20260914T193517Z/tables/c1_contrasts.tsv",
        PROG / "c1-matched-endogenous-20260914T193517Z/tables/c1_model_metrics.tsv",
        PROG / "c1-matched-endogenous-20260914T193517Z/tables/c1_matched_counts.tsv",
        PROG / "c1-endpoint2-gpu-20260914T203750Z/tables/c1e2_contrasts.tsv",
        PROG / "c1-endpoint2-gpu-20260914T203750Z/tables/c1e2_model_metrics.tsv",
        PROG / "c1-endpoint2-gpu-20260914T203750Z/tables/c1e2_matched_counts.tsv",
        PROG / "c2-endogenous-head-20260914T194000Z/tables/head_summary.tsv",
        PROG / "c2-endogenous-head-20260914T194000Z/tables/learning_curve.tsv",
        PROG / "d1-lineage-allelic-20260914T194213Z/tables/lineage_feasibility.tsv",
        PROG / "e1-bridge-20260914T194145Z/bridge_counts.tsv",
        PROG / "e1-bridge-20260914T194145Z/e2_distance_strata.tsv",
        PROG / "e2-transfer-20260914T203752Z/tables/defined_contrast_summary.tsv",
        PROG / "e2-transfer-20260914T203752Z/tables/measured_label_ceiling.tsv",
        PROG / "e2-transfer-20260914T203752Z/tables/positive_control.tsv",
        PROG / "e3-reporter-context-20260914T193515Z/tables/e3_correlations.tsv",
        PROG / "e3-reporter-context-20260914T193515Z/tables/e3_predictions.tsv",
        PROG / "f1-phase-20260914T193608Z/tables/f1_summary.json",
        PROG / "f1-phase-20260914T193608Z/tables/f1_read_backed_phase_summary.tsv",
        PROG / "f1-phase-20260914T193608Z/tables/topld_allele_concordance.tsv",
        PROG / "f2-haplotype-v2-20260915T063000Z/tables/f2v2_summary.json",
        PROG / "f2-haplotype-v2-20260915T063000Z/tables/f2v2_gnmt_verdict.json",
        PROG / "b1-independent-check-20260914T204500Z/check_b1_results.json",
        PROG / "b1-v2-independent-check-20260914T214500Z/check_b1_v2_results.json",
        PROG / "e3-independent-check-20260914T213000Z/check_results.json",
        ROOT / "scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md",
    ]
    rows = []
    for p in inputs:
        rows.append(
            {
                "path": str(p),
                "exists": p.exists(),
                "bytes": p.stat().st_size if p.exists() else "",
                "sha256": sha256(p) if p.exists() else "",
            }
        )
    rows.append(
        {
            "path": str(Path(__file__).resolve()),
            "exists": True,
            "bytes": Path(__file__).stat().st_size,
            "sha256": sha256(Path(__file__)),
        }
    )
    pd.DataFrame(rows).to_csv(OUT / "MANIFEST.tsv", sep="\t", index=False)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
