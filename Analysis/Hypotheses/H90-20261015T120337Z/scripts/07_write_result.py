"""H90 final step when a stop rule fires before the outcome: compute the verdict from the gate
files and write RESULT_H90.md. Reads no program score (none exists for GSE289173)."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import h90_common as C

C.check_inputs()
out = C.OUT
res = C.RUN / "RESULT_H90.md"
if res.exists():
    sys.exit("refusing to overwrite RESULT_H90.md")
if (out / "donor_scores.tsv").exists() or (out / "results_primary.tsv").exists():
    sys.exit("an outcome file exists; this script is only for a pre-outcome stop")

g0 = pd.read_csv(out / "gate_G0_download.tsv", sep="\t")
g4 = pd.read_csv(out / "gate_G4_bridge.tsv", sep="\t")
g13 = pd.read_csv(out / "gate_G1_G3.tsv", sep="\t", index_col=0)["value"]
rc = pd.read_csv(out / "G3_recheck_summary.tsv", sep="\t", index_col=0)
genes = pd.read_csv(out / "G3_recheck_genes.tsv", sep="\t")
d4 = pd.read_csv(out / "D4_bridge.tsv", sep="\t")
pc2 = pd.read_csv(out / "PC2_collagen_result.tsv", sep="\t").iloc[0]
tobool = lambda v: str(v).strip().lower() == "true"

g0_pass = g0["pass"].map(tobool).all()
g4a = g4[g4.check.str.startswith("G4a")]["pass"].map(tobool).all()
g4b = tobool(g4.loc[g4.check == "G4b_projection_slope_positive", "pass"].iloc[0])
g1 = tobool(g13["G1_pass"]); g2 = tobool(g13["G2_pass"])
prim_scored = int(float(g13["primary_scored"])); prim_w = float(g13["primary_weight"])
g3 = prim_scored >= 8 and prim_w >= 0.80
# the independent recheck must agree with job 6301
rc20 = rc.loc[20]
if int(rc20.n_scored) != prim_scored or abs(float(rc20.retained_weight) - prim_w) > 1e-9:
    sys.exit("G3 recheck disagrees with job 6301; resolve before writing a verdict")

order = [("G0", g0_pass, "untestable: deposited matrices do not match the metadata"),
         ("G1", g1, "untestable: sample labels cannot be linked to matrices"),
         ("G2", g2, "untestable: too few donors"),
         ("G3", g3, "untestable: program not measured"),
         ("G4a", g4a, "untestable: Fig 4C universe or refit not reproduced"),
         ("G4b", g4b, "untestable as a replication: the frozen-weight projection does not carry the Fig 4C direction in discovery")]
failed = [(k, msg) for k, ok, msg in order if not ok]
if not failed:
    sys.exit("no gate failed; run 04_score_programs.py and 05_fit_h90.R instead")
verdict = failed[0][1]
sec_scored = int(float(g13["secondary_scored"])); sec_w = float(g13["secondary_weight"])
low = genes[(genes.module == 20) & ~genes.scored].sort_values("detection")
lowtxt = ", ".join(f"{g} {100 * d:.2f}%" for g, d in zip(low.gene, low.detection))
dp = d4[(d4.score == "projection") & (d4.dataset == "n-weighted")].iloc[0]
df4 = d4[(d4.score == "fig4c_hotspot") & (d4.dataset == "n-weighted")].iloc[0]
rho = float(g4.loc[g4.check.str.startswith("G4c"), "value"].iloc[0])

md = f"""# RESULT_H90

Run `{C.RUN.name}`. Prespec sha256 {C.PRESPEC_SHA}. Verdict computed in `scripts/07_write_result.py` from the gate files.

**Primary verdict (Ductular injury, hepatocyte module 20): {verdict.upper()}** (first failed gate: {', '.join(k for k, _ in failed)})

No GSE289173 program score was computed and no outcome test was run. Under the prespec "Untestable" reading, no H90 number is used in any claim.

Supplement sentence (the only text the prespec allows): GSE289173 (48 donors, single-nucleus RNA-seq) was opened to test whether the Ductular injury score is higher in later-stage donors, but the test could not be run because only {prim_scored} of 21 program genes, carrying {prim_w:.2f} of the frozen program weight, were detected in at least 1% of hepatocyte nuclei (prespecified minimum 0.80).

## Gates

| Gate | Result |
|---|---|
| G0 data identity | {"pass" if g0_pass else "FAIL"} (48 matrices; 283,890 barcodes match the metadata; features identical) |
| G1 labels (PC1) | {"pass" if g1 else "FAIL"}; {int(float(g13["n_discordant_sex"]))} discordant donor dropped ({g13["discordant"]}) |
| G2 minimum n | {"pass" if g2 else "FAIL"}; ctrl {int(float(g13["n_ctrl"]))}, MASL {int(float(g13["n_MASL"]))}, MASH {int(float(g13["n_MASH"]))} |
| G3 coverage, Ductular injury | {"pass" if g3 else "FAIL"}; {prim_scored} genes, retained weight {prim_w:.3f} (needs >= 8 and >= 0.80); independently re-derived from the raw GEO matrices |
| G3 coverage, Stromal ECM | {sec_scored} genes, weight {sec_w:.3f}; secondary untestable |
| G4 bridge | {"pass" if (g4a and g4b) else "FAIL"}; Fig 4C refit reproduced; projection slope in discovery > 0; within-dataset Spearman with the Fig 4C score {rho:.2f} |
| G5 | not reached |

Ductular injury genes below 1% detection in GSE289173 hepatocyte nuclei: {lowtxt}.

## Pre-outcome diagnostics (description only, not for claims)

- Discovery bridge (D4): standardized stage means of the frozen-weight projection {dp.m1:.2f} (Steatosis) and {dp.m2:.2f} (Steatohepatitis), against {df4.m1:.2f} and {df4.m2:.2f} for the Fig 4C Hotspot score; implied slope at the planned 23/6/19 design {dp.implied_slope_planned_23_6_19:.2f} SD per step, expected power {dp.expected_power_alpha05_planned:.2f} at alpha 0.05 (Fig 4C score {df4.implied_slope_planned_23_6_19:.2f}, power {df4.expected_power_alpha05_planned:.2f}).
- PC2 (collagen in fibroblast nuclei, aMASH vs MASL): +{pc2.beta_aMASH_vs_MASL:.2f} log2 CPM, standardized difference {pc2.std_diff:.2f} ({int(pc2.n_MASL)} MASL, {int(pc2.n_aMASH)} aMASH donors).
- Realized design: r(stage, log10 hepatocyte nuclei) {float(g13["r_stage_log10nhep"]):.2f}; VIF {float(g13["VIF_stage"]):.2f}; MDE {float(g13["MDE_realized_x_sqrtVIF"]):.2f} SD per step.

## Reading under the H59 table

Use the "H90 untestable" column of the PRESPEC_H90 "Dependence on H59" table. Stromal ECM also takes the "H90 untestable" column.

## Sealed status

GSE289173 is registered `withheld_sealed` in `Analysis/MASLD_Model_Benchmark/config/datasets/gse289173.toml`. Two events ended that status: on 2026-10-15 the feasibility gate read disease group and sex for all 48 samples from the GEO series matrix (beyond the 2 samples in the benchmark's 2026-08-21 incident), and on 2026-10-15 this run downloaded and opened all 48 expression matrices (QC, CellTypist labels, sex genes, collagen and program-gene detection). This run did not edit the registry; the owner should update it.
"""
res.write_text(md)
print(md)
