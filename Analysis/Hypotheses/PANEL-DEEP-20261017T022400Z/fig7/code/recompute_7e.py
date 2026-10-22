"""Recompute every Figure 7E row (observed / its own decision threshold) from the upstream result files,
with the definition that reproduces the drawn value and the alternatives a reader could reasonably expect.
Uses this review's own eQTL strata and routing recomputations for those two rows."""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
EX = ROOT / "Analysis/MASLD_Model_Benchmark/executions"
OUT = Path(sys.argv[1])
rows = []


def add(label, verdict_drawn, value, definition, alt=None, alt_def=None, alt_verdict=None, n=None):
    rows.append({"row": label, "verdict_drawn": verdict_drawn, "ratio": value, "definition": definition,
                 "alternative": alt, "alternative_definition": alt_def, "verdict_under_alternative": alt_verdict, "n": n})


arms = {a["arm"]: a for a in json.load(open(EX / "crossmodal-severity-stratification-20260831T190048Z/results/ALL_ARMS_SUMMARY.json"))}
p = arms["primary_activity_ge3"]
add("normal vs diseased map", "supported", p["delta"] / p["mde"], "GSE267145 primary_activity_ge3 delta / MDE",
    p["centred"] / p["mde"], "null-centred delta / MDE", "supported", f"{p['n_low']}+{p['n_high']} participants")
z = json.load(open(EX / "gse105127-zonation-lane-20260830T195218Z/results/zone_metrics.json"))
acc = z["primary_result_free_argmax"]["accuracy"]; bars = z["preregistered_bars"]
add("zone from expression", "supported", acc / bars["significance_bar_accuracy_0.05"], "accuracy / significance-bar accuracy (0.4737)",
    (acc - bars["chance_accuracy"]) / (bars["minimum_detectable_effect_accuracy"] - bars["chance_accuracy"]),
    "(accuracy - chance) / (MDE accuracy - chance); chance alone scores 0.70 on the drawn scale", "supported", "19 participants, 57 rows")
tr = pd.read_csv(ROOT / "GWAS/ldsc/program_heritability_20260830/config/traits.tsv", sep="\t", comment="#").set_index("trait")
H2BAR = 7.0
add("NAFLD endpoint", "untestable", tr.loc["finngen_nafld", "h2_z"] / H2BAR, "FinnGen R12 NAFLD h2 z / 7",
    tr.loc["mvp_nafld", "h2_z"] / H2BAR, "MVP NAFLD h2 z / 7 (largest NAFLD GWAS h2 z)", "untestable")
b = arms["posthoc_B_female_only_NOR_dropped_activity_ge4"]; a = arms["posthoc_A_female_only_activity_ge3"]
add("severity, female only", "untestable", b["delta"] / b["mde"], "GSE267145 post hoc B: female only AND normal livers dropped AND activity >= 4",
    a["delta"] / a["mde"], "post hoc A: female only at the primary threshold (activity >= 3)", a["verdict"], f"B {b['n_low']}+{b['n_high']}; A {a['n_low']}+{a['n_high']}")
f = arms["primary_fibrosis_any"]
add("ATAC map vs fibrosis", "untestable", f["delta"] / f["mde"], "GSE296875 primary_fibrosis_any delta / MDE", None, None, None, f"{f['n_low']}+{f['n_high']} donors")
add("NAFLD meta-analysis", "untestable", tr.loc["ghodsian_nafld", "h2_z"] / H2BAR, "Ghodsian 2021 h2 z / 7")
w = arms["posthoc_NOR_dropped_activity_ge4"]
add("severity within disease", "untestable", w["delta"] / w["mde"], "GSE267145 post hoc: normal livers dropped, activity >= 4", None, None, None, f"{w['n_low']}+{w['n_high']}")
add("cirrhosis endpoint", "untestable", tr.loc["finngen_cirr", "h2_z"] / H2BAR, "FinnGen R12 cirrhosis h2 z / 7",
    tr.loc["ghouse_cirr", "h2_z"] / H2BAR, "Ghouse 2023 cirrhosis h2 z / 7 (run as secondary trait, 0 of 117 programs at BH 0.05)", "untestable")
s = arms["secondary_steatosis_ge5pct"]
add("ATAC map vs steatosis", "untestable", s["delta"] / s["mde"], "GSE296875 secondary_steatosis_ge5pct delta / MDE", None, None, None, f"{s['n_low']}+{s['n_high']}")
st = pd.read_csv(OUT / "eqtl_strata_recomputed.tsv", sep="\t").set_index("stratum").loc["all"]
add("signed eQTL direction", "tested negative", (st.best_auc - st.chance_median) / st.resolution_half_width,
    "(best AUC - within-block chance) / CI half-width, all rows (this review's bootstrap)",
    st.seq_bound_minus_chance / (2 * st.resolution_half_width),
    "(sequence-attributable bound - chance) / (2 x resolution), the panel's own untestable rule applied to the bound drawn in 7D",
    "untestable", f"{int(st.heldout_rows)} rows, {int(st.heldout_blocks)} LD blocks")
r = json.load(open(OUT / "fibrosis_routing_rrr_recomputed.json"))["routing"]
add("routing adds skill", "tested negative", r["oracle_gain_over_floor_ratio"], "oracle-router gain over reference / prespecified floor 0.093",
    None, None, None, "521 samples, 5 cohorts, no donor key")
ph = pd.read_csv(ROOT / "GWAS/ldsc/program_heritability_20260830/results/program_heritability_summary.tsv", sep="\t")
for trait, lab in [("ukbb_ast", "program h2, AST"), ("ukbb_alt", "program h2, ALT")]:
    g = ph[ph.trait == trait]
    add(lab, "tested negative", float((g.tau_p_bh_within_trait < 0.05).sum()), "programs with conditional tau BH < 0.05 (of 117)",
        float((g.enrichment_p_bh_within_trait < 0.05).sum()), "programs with marginal enrichment BH < 0.05 (of 117)", "not a verdict change; different test", "117 programs")
wa = json.load(open(EX / "gse83452-intervention-response-20260831/results_20260831/within_arm_decisive.json"))["readouts"]["gene_median_elastic_net"]["diet_arm_only"]
add("treatment response", "tested negative", -wa["diff"] / wa["mde_80pct"], "-(diet-arm resolved minus persisted) / MDE80, oriented to prespecified direction",
    None, None, None, f"{wa['n_resolved']} resolved vs {wa['n_persisted']} persisted")
t = pd.DataFrame(rows)
t.to_csv(OUT / "panel7e_recomputed.tsv", sep="\t", index=False, float_format="%.4f")
print(t[["row", "verdict_drawn", "ratio", "alternative", "verdict_under_alternative"]].to_string())
