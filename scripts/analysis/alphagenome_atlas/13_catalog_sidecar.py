#!/usr/bin/env python3
"""Step 13: Gene Catalog sidecar keyed on ensembl_id + credible_set_id (v2 contract untouched).

Fields: candidate mechanism (channel families with strong liver prediction; never a single argmax),
measured support (chromatin/protein/remodeling states), unresolved alternative, next experiment rule id
(reuses the v2 rule ids; deterministic, outcome-blind routing).
"""

from __future__ import annotations

import math

import pandas as pd

import lib_atlas as la

TABLES = la.out_root() / "tables"
STRONG = 0.9


def main() -> None:
    prof = pd.read_csv(TABLES / "signal_profiles.tsv", sep="\t", keep_default_na=False)
    chan = pd.read_csv(TABLES / "signal_channel_profiles.tsv", sep="\t", keep_default_na=False)
    census = pd.read_csv(TABLES / "chromatin_census.tsv", sep="\t", keep_default_na=False).set_index("signal_uid")
    ctx = pd.read_csv(TABLES / "context_summary.tsv", sep="\t", keep_default_na=False).set_index("signal_uid") if (TABLES / "context_summary.tsv").exists() else None
    dirs = pd.read_csv(TABLES / "direction_comparisons.tsv", sep="\t", keep_default_na=False)
    d1 = dirs[dirs["comparison"] == "predicted_vs_measured_eqtl"].set_index("signal_uid")
    rows = []
    for _, p in prof.iterrows():
        sig = p["signal_uid"]
        c = chan[(chan["signal_uid"] == sig) & (chan["liver_tracks_available"].astype(str) == "True")]
        strong = [r["channel"] for _, r in c.iterrows() if _f(r["A_quantile_over_C"]) >= STRONG and _f(r["C_liver"]) >= 0.5]
        weak_cov = [r["channel"] for _, r in c.iterrows() if _f(r["C_liver"]) < 0.5]
        if p["universe"] == "B_direct":
            cs_id = sig.replace("susie:", "")
        else:
            cs_id = sig.replace("coloc:", "coloc:")
        cc = census.loc[sig] if sig in census.index else None
        measured = []
        if cc is not None:
            for lin in ("hepatocyte", "stellate"):
                for st in ("supported", "discordant", "source_dependent"):
                    m = _f(cc.get(f"{lin}_da_{st}_mass", math.nan))
                    if m >= 0.5:
                        measured.append(f"{lin}_accessibility_change_{st}(mass={m:.2f})")
        if ctx is not None and sig in ctx.index:
            x = ctx.loc[sig]
            if x["remodeling_state"]:
                measured.append(f"bulk_remodeling={x['remodeling_state']}")
            if x["protein_state"]:
                measured.append(f"protein={x['protein_state']}")
        coding = _f(p["coding_mass"])
        mech = []
        if coding >= 0.5:
            mech.append("protein_altering")
        mech += [f"predicted_{s}" for s in strong]
        if not mech:
            mech.append("unresolved_noncoding" if coding < 0.5 else "protein_altering")
        alt = []
        if str(p["flag_multiple_candidate_genes"]) == "True":
            alt.append(f"other_gene_expression_prediction:{p['strong_expression_genes']}")
        if str(p["flag_opposing_signs_target_expression"]) == "True":
            alt.append("credible_variants_predict_opposite_expression_directions")
        if str(p["flag_quiet_expression_with_processing"]) == "True":
            alt.append("processing_rather_than_abundance")
        if str(p["flag_strong_avi_little_liver"]) == "True":
            alt.append("generic_impact_without_liver_prediction")
        if weak_cov:
            alt.append("low_liver_coverage:" + ";".join(weak_cov))
        if sig in d1.index and d1.at[sig, "state"] in ("unresolved", "discordant"):
            alt.append(f"eqtl_direction_{d1.at[sig, 'state']}")
        if coding >= 0.5:
            rule = "EXP_PROTEIN_STATE_CONTEXT_V2"
        elif strong:
            rule = "EXP_REGULATORY_DNA_V2"
        elif not c.empty and all(_f(r["C_liver"]) < 0.5 for _, r in c.iterrows()):
            rule = "EXP_MEASURE_MISSING_ASSAY_V2"
        else:
            rule = "EXP_REGULATORY_DNA_V2"
        rows.append({"ensembl_id": p["ensembl"], "credible_set_id": cs_id, "signal_uid": sig, "universe": p["universe"], "gwas_name": p["gwas_name"],
                     "trait_class": p["trait_class"], "analysis_block": p["analysis_block"], "posterior_definition": p["posterior_definition"],
                     "candidate_mechanism": ";".join(mech), "coding_mass": coding, "consequence_coverage": p["consequence_C"],
                     "avi_top_quantile": p["avi_top_quantile"], "avi_top_features": p["avi_top_features"],
                     "measured_support": ";".join(measured) if measured else "none_recorded",
                     "unresolved_alternative": ";".join(alt) if alt else "none_flagged", "recommended_experiment_rule_id": rule, "next_experiment_rule_id": rule,   # v2 field name plus the legacy alias
                     "atlas_extension_state": "candidate_not_adopted"})
    la.write_tsv_once(TABLES / "gene_catalog_atlas_sidecar.tsv", rows, list(rows[0].keys()))
    la.log(f"sidecar: {len(rows)} rows")


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


if __name__ == "__main__":
    main()
