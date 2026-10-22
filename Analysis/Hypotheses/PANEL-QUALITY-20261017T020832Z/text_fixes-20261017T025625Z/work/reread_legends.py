"""Re-read every number used in the L- edits to docs/manuscript/05_figure_legends.md
from its producing file. Read-only on all inputs.
Writes TSV to stdout: edit_id, source_path, column_or_row, value_reread."""
import hashlib
import json
import math
import os
import re
import subprocess

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
rows = []


def p(rel):
    return os.path.join(ROOT, rel)


def add(eid, path, col, val):
    rel = os.path.relpath(path, ROOT) if path.startswith(ROOT) else path
    rows.append((eid, rel, col, str(val)))


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def grep(path, pattern):
    with open(path) as fh:
        return [ln.strip() for ln in fh if re.search(pattern, ln)]


def rds_dims(path):
    cmd = [RSCRIPT, "-e", f'suppressMessages(library(edgeR)); d <- readRDS("{path}"); cat(dim(d$counts))']
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=300).stdout.split()
    return tuple(int(x) for x in out[-2:])


# ---- L-01 Fig 1B: 774-unit interface is genetic-identity clusters
arm = p("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/arm_deg_counts.tsv")
a = pd.read_csv(arm, sep="\t").set_index("arm").loc["primary_source_controls"]
for c in ("n_samples", "n_control", "n_disease", "n_genes"):
    add("L-01", arm, f"primary_source_controls {c}", int(a[c]))
ros = p("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_unit_roster.tsv")
r = pd.read_csv(ros, sep="\t")
add("L-01", ros, "n rows; n unique analysis_unit_id", f"{len(r)}; {r['analysis_unit_id'].nunique()}")
add("L-01", ros, "unique unit_evidence", "|".join(sorted(r["unit_evidence"].astype(str).unique())))

# ---- L-02 Fig 2A: PP.H4 bands
casc = p("figures/main/fig2_genetics/panels/fig2A_gwas_cascade_source.tsv")
c = pd.read_csv(casc, sep="\t")
pb = c[(c["record_type"] == "posterior_band") & (c["metric"] == "coloc_genes")]
for _, x in pb.iterrows():
    add("L-02", casc, f"posterior_band {x['signal_model']} {x['pp_h4_band']} value", int(x["value"]))
for sm, d in pb.groupby("signal_model"):
    add("L-02", casc, f"posterior_band {sm} sum(value)", int(d["value"].sum()))
gw = p("scripts/figures/fig1_gwas_creative_options.py")
for ln in grep(gw, r'su_ok, abf_ok = |band = ">0\.9"|> g_bs\[gene\]|> g_ba\[gene\]'):
    add("L-02", gw, "band rule line", ln)

# ---- L-03 Fig 2B: dashed line and labels
inst2b = p("figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc.pdf")
rep2b = p("figures/candidates/review-r3-fixes-20261015T112148Z/F2abce/repair/Fig2B_PIP_vs_SuSiE-coloc_sameRow_susie_repair.pdf")
add("L-03", inst2b, "sha256 installed == repair PDF", sha(inst2b) == sha(rep2b))
src2b = p("figures/candidates/review-r3-fixes-20261015T112148Z/F2abce/repair/Fig2B_PIP_vs_SuSiE-coloc_sameRow_susie_repair_source.csv")
b = pd.read_csv(src2b)
add("L-03", src2b, "n rows (plotted genes)", len(b))
add("L-03", src2b, "genes with is_label TRUE", "|".join(sorted(b.loc[b["is_label"], "gene"])))
roster = ["HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "RORA"]
add("L-03", src2b, "roster genes present among plotted rows", "|".join(g for g in roster if g in set(b["gene"])))
r2b = p("Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce/scripts/fig2b_render_repair.R")
old2b = p("figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv")
ob = pd.read_csv(old2b)
add("L-03", old2b, "installed source CSV: n rows; n is_label (does not match installed PDF)", f"{len(ob)}; {int(ob['is_label'].sum())}")
for ln in grep(r2b, r"^roster <-|geom_hline|aes\(x = x_max, y = pp4p\)"):
    add("L-03", r2b, "render line", ln)

# ---- L-04 Fig 2G-H: PIP track source
fz = p("scripts/figures/figS09_locus_zoom_disease.R")
for f in ("Fig2G_RORA_NAFLD_source.csv", "Fig2H_FABP1_ALT_source.csv"):
    path = p("figures/main/fig2_genetics/panels/" + f)
    g = pd.read_csv(path).iloc[0]
    add("L-04", path, "gene; gwas; gwas_ancestry; pip_source; n_pip_variants; plot_ld_source",
        f"{g['gene']}; {g['gwas']}; {g['gwas_ancestry']}; {g['pip_source']}; {int(g['n_pip_variants'])}; {os.path.basename(g['plot_ld_source'])}")
for ln in grep(fz, r'pip_source <- "|susie_rss\(z = z, R = R, L = 1, coverage = 0\.95'):
    add("L-04", fz, "PIP source line", ln)

# ---- L-05 Fig 3G: lncRNA conditional OR and interval type
s3 = p("figures/main/fig3_bulk_transcriptomics/source_tables/fig3_lncrna_role/summary.tsv")
s = pd.read_csv(s3, sep="\t")
md = s[(s["section"] == "matched_direction")].set_index("metric")
orr = md.loc["conditional_odds_ratio"]
add("L-05", s3, "matched_direction conditional_odds_ratio value; lower; upper; n",
    f"{float(orr['value']):.6f}; {float(orr['lower']):.6f}; {float(orr['upper']):.6f}; {int(float(orr['n']))}")
nl = int(float(md.loc["discordant_lncrna_up", "value"]))
nc = int(float(md.loc["discordant_coding_up", "value"]))
add("L-05", s3, "discordant_lncrna_up; discordant_coding_up", f"{nl}; {nc}")
se = math.sqrt(1 / nl + 1 / nc)
lo, hi = (math.exp(math.log(nl / nc) + k * 1.96 * se) for k in (-1, 1))
add("L-05", s3, "recomputed OR = b/c; Wald CI exp(log OR +/- 1.96*sqrt(1/b+1/c))", f"{nl / nc:.6f}; {lo:.6f}; {hi:.6f}")
r3 = p("scripts/figures/fig3_lncrna_role.R")
for ln in grep(r3, r"matched_ci <- exp"):
    add("L-05", r3, "CI line", ln)

# ---- L-06 Fig 4E: covariates and HC3
inst4e = p("figures/main/fig4_singlecell_programs/panels/fig4e_bulk_tissue_state_transport.pdf")
cand4e = p("figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4e/panels/fig4e_bulk_tissue_state_transport.pdf")
add("L-06", inst4e, "sha256 installed == review-r3 fig4e PDF", sha(inst4e) == sha(cand4e))
pe = p("figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4e/source_tables/fig4e_plotted_effects.tsv")
add("L-06", pe, "unique n_cohorts", "|".join(str(x) for x in sorted(pd.read_csv(pe, sep="\t")["n_cohorts"].unique())))
ce = p("figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_three_cohorts/cohort_effects.tsv")
e = pd.read_csv(ce, sep="\t")
prim = e[e["model_kind"] == "primary"]
add("L-06", ce, "se_hc3 column present", "se_hc3" in e.columns)
add("L-06", ce, "primary unique n_model_columns", "|".join(str(int(x)) for x in sorted(prim["n_model_columns"].dropna().unique())))
add("L-06", ce, "primary cohorts", "|".join(sorted(prim["cohort"].unique())))
pr = p("figures/candidates/fig4e-program-axis-774-20261012T010830Z/producers/program_axis_map_774.R")
for ln in grep(pr, r'common <- c\("sample_id"'):
    add("L-06", pr, "model columns line", ln)

# ---- L-07 Fig 4F: windows and per-cohort n
inst4f = p("figures/main/fig4_singlecell_programs/panels/fig4f_continuum_program_trajectories.pdf")
cand4f = p("figures/candidates/fig3-fig4-layout-774-20261013T075956Z/panels/fig4/fig4f_continuum_program_trajectories.pdf")
if os.path.exists(inst4f):
    add("L-07", inst4f, "sha256 installed == layout-774 candidate PDF", sha(inst4f) == sha(cand4f))
src4f = p("figures/candidates/fig3-fig4-layout-774-20261013T075956Z/panels/fig4/fig4f_continuum_program_trajectories_source.tsv")
w = pd.read_csv(src4f, sep="\t")
win = w.drop_duplicates("window_id").sort_values("window_id")
add("L-07", src4f, "n windows", len(win))
add("L-07", src4f, "window_lower..window_upper", "|".join(f"{a:.1f}-{b:.1f}" for a, b in zip(win["window_lower"], win["window_upper"])))
for ds, d in w.groupby("dataset"):
    add("L-07", src4f, f"{ds} n_participants min-max", f"{d['n_participants'].min()}-{d['n_participants'].max()}")
ps = p("RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261011T152106Z/projection/participant_scores.tsv")
q = pd.read_csv(ps, sep="\t", usecols=["sample_id", "dataset"])
n_by = q.groupby("dataset")["sample_id"].nunique()
for ds, n in n_by.items():
    add("L-07", ps, f"{ds} unique sample_id", int(n))
add("L-07", ps, "total unique sample_id", int(q["sample_id"].nunique()))
fp = pd.read_csv(ps, sep="\t", usecols=["dataset", "axis_id", "axis_percentile"])
fp = fp[fp["axis_id"] == "fixed_projection"]
for ds, d in fp.groupby("dataset"):
    add("L-07", ps, f"{ds} fixed_projection axis_percentile min-max (within cohort)", f"{d['axis_percentile'].min():.2f}-{d['axis_percentile'].max():.2f}")

# ---- L-08 Fig 4G: expression universe of scores vs geometry
ml = p("RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/inputs/frozen_contract.json")
add("L-08", ml, "gene_family_size", json.load(open(ml))["gene_family_size"])
fh = p("Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F4efgh/fig4h_disease_only/fig4h_disease_only.py")
for ln in grep(fh, r"^ML = "):
    add("L-08", fh, "score input line", ln)
dge774 = p("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds")
add("L-08", dge774, "dim(counts) genes x units", "x".join(map(str, rds_dims(dge774))))
gim = p("figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/provenance/geometry_input_manifest.tsv")
gm = pd.read_csv(gim, sep="\t")
eu = gm.loc[gm["input_role"] == "expression_universe_only", "path"].iloc[0]
add("L-08", gim, "expression_universe_only path", eu.replace("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/", ""))
add("L-08", eu, "dim(counts) genes x records", "x".join(map(str, rds_dims(eu))))
gp = p("figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/geometry/geometry_parameters.tsv")
g = pd.read_csv(gp, sep="\t").set_index("parameter")["value"]
add("L-08", gp, "seed; gene_count; feature_count", f"{g['seed']}; {g['gene_count']}; {g['feature_count']}")
inst4g = p("figures/main/fig4_singlecell_programs/panels/fig4g_molecular_systems_continuum_effect_map.pdf")
cand4g = p("figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/repair/fig4g/fig4g_molecular_systems_continuum_effect_map_disease_only.pdf")
add("L-08", inst4g, "sha256 installed == review-r3 disease-only 4G", sha(inst4g) == sha(cand4g))

# ---- L-09 Fig 5B: protein significance rule and BH family
pde = p("Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv")
d = pd.read_csv(pde, sep="\t")
bh = multipletests(d["pvalue"], method="fdr_bh")[1]
add("L-09", pde, "n proteins (rows); unique gene", f"{len(d)}; {d['gene'].nunique()}")
add("L-09", pde, "max |BH(pvalue) - padj|", f"{np.abs(bh - d['padj']).max():.2e}")
f5b = p("scripts/figures/fig5b_protein_triage.R")
for ln in grep(f5b, r"protein_padj < 0\.05, \"confirmed\""):
    add("L-09", f5b, "significance rule line", ln)
sm = p("figures/main/fig5_molecular_context/panels/data/fig5b_protein_triage_summary.tsv")
sv = pd.read_csv(sm, sep="\t").iloc[0]
add("L-09", sm, "n_protein_significant; n_prioritized_measured", f"{sv['n_protein_significant']}; {sv['n_prioritized_measured']}")

# ---- L-10 Fig 5C: histology BH family and enrichment test
B = "figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/"
hp = p(B + "panel4c_histology_partial.tsv")
h = pd.read_csv(hp, sep="\t")
add("L-10", hp, "n rows; n genes; features", f"{len(h)}; {h['gene'].nunique()}; {'|'.join(sorted(h['feature'].unique()))}")
add("L-10", hp, "unique correction_family", "|".join(h["correction_family"].unique()))
add("L-10", hp, "max |BH(pvalue) - padj| over all rows", f"{np.abs(multipletests(h['pvalue'], method='fdr_bh')[1] - h['padj']).max():.2e}")
en = p(B + "protein_pathway_enrichment.tsv")
x = pd.read_csv(en, sep="\t")
add("L-10", en, "n sets; family counts", f"{len(x)}; " + "; ".join(f"{k} {v}" for k, v in x["family"].value_counts().items()))
add("L-10", en, "max |BH(pval) - padj|", f"{np.abs(multipletests(x['pval'], method='fdr_bh')[1] - x['padj']).max():.2e}")
pp = p("Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R")
for ln in grep(pp, r"ranks <- sort\(setNames\(ranked\$t|minSize = 8|maxSize = 600|^set\.seed"):
    add("L-10", pp, "fgsea line", ln)
pc = p(B + "panel4c_process_enrichment.tsv")
z = pd.read_csv(pc, sep="\t")
for _, y in z[z["representable"]].sort_values("display_order").iterrows():
    add("L-10", pc, f"{y['program']} NES; padj", f"{y['NES']:.2f}; {y['padj']:.1e}")
    assert y["pathway"] in set(x["pathway"]), y["pathway"]

# ---- L-11 Fig 5D: bootstrap draws and seed
s5 = p("figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_locus_background_summary.tsv")
t = pd.read_csv(s5, sep="\t")
add("L-11", s5, "unique bootstrap", "|".join(t["bootstrap"].unique()))
le = t[t["trait_class"] == "liver_enzyme"]
add("L-11", s5, "liver_enzyme n_genes; min ratio_lo", f"{le['n_genes'].iloc[0]}; {le['ratio_lo'].min():.4f}")
f5d = p("figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/scripts/fig5d_summaries.py")
for ln in grep(f5d, r"default_rng|rng\.integers|ratio_lo=np\.percentile"):
    add("L-11", f5d, "bootstrap line", ln)

# ---- L-12 Fig 7B: held-out fold, n, paired bootstrap
iv = p("Analysis/MASLD_Model_Benchmark/executions/dnalm_heads-20260831T183647Z/logs/03_intervals.json")
j = json.load(open(iv))
add("L-12", iv, "seed; bootstrap_resamples; bootstrap_unit; primary_held_out_fold",
    f"{j['seed']}; {j['bootstrap_resamples']}; {j['bootstrap_unit']}; {j['primary_held_out_fold']}")
g5 = j["G5_baseline_reproduction"]
add("L-12", iv, "heldout_rows_here; heldout_blocks_here", f"{g5['heldout_rows_here']}; {g5['heldout_blocks_here']}")
for m in ("hyenadna", "nucleotide_transformer"):
    k = f"{m}::m1_dist_emb  MINUS  {m}::b1_dist"
    dd = j["paired_deltas"]["abs_z"][k]
    add("L-12", iv, f"paired_deltas abs_z {k} point; ci", f"{dd['point']:.6f}; {dd['ci'][0]:.6f} to {dd['ci'][1]:.6f}")
g7b = p("Analysis/MASLD_Model_Benchmark/executions/fig7_panels_rebuild-20260901T012245Z/tables/panel_7b_gains.tsv")
gb = pd.read_csv(g7b, sep="\t")
gb = gb[gb["task"].str.startswith("variant -> gene |z|")]
add("L-12", g7b, "task_metric; unit_and_n; interval_kind",
    f"{gb['task_metric'].iloc[0]}; {gb['unit_and_n'].iloc[0]}; {gb['interval_kind'].iloc[0]}")
s3i = p("Analysis/MASLD_Model_Benchmark/executions/dnalm_heads-20260831T183647Z/scripts/03_intervals.py")
for ln in grep(s3i, r"np\.quantile"):
    add("L-12", s3i, "percentile line", ln)

# ---- L-13 Fig 7C: no Geneformer row in the installed panel
inst7c = p("figures/main/fig7_model_framework/fig7C_score_beside_training_exposure.pdf")
cand7c = p("Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/panels/fig7C_score_beside_training_exposure.pdf")
add("L-13", inst7c, "sha256 installed == exposure-correction 7C", sha(inst7c) == sha(cand7c))
t7c = p("Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/tables/panel_7c_exposure.tsv")
c7 = pd.read_csv(t7c, sep="\t")
add("L-13", t7c, "rows whose model contains 'geneformer'", int(c7["model"].str.contains("geneformer", case=False).sum()))
txt7c = subprocess.run(["pdftotext", inst7c, "-"], capture_output=True, text=True).stdout
add("L-13", inst7c, "pdftotext contains 'Geneformer'", "geneformer" in txt7c.lower())

# ---- L-14 / L-15 Fig 7D: interval method and block n
inst7d = p("figures/main/fig7_model_framework/fig7D_measured_bounds_and_join_structure.pdf")
cand7d = p("figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair/fig7D_measured_bounds_and_join_structure.pdf")
add("L-14", inst7d, "sha256 installed == F7b repair 7D", sha(inst7d) == sha(cand7d))
tu = p("figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair/fig7_transfer_primary_units.tsv")
u = pd.read_csv(tu, sep="\t")
add("L-14", tu, "unique ci_B; ci_seed; ci_method",
    f"{'|'.join(map(str, u['ci_B'].unique()))}; {'|'.join(map(str, u['ci_seed'].unique()))}; {'|'.join(u['ci_method'].unique())}")
add("L-15", iv, "signed-eQTL held-out rows; LD blocks", f"{g5['heldout_rows_here']}; {g5['heldout_blocks_here']}")
E = p("Analysis/MASLD_Model_Benchmark/executions/variant_program_substrate_rebuild-20260830T160000Z/")
tm = os.path.join(E, "tables/training_matrix.parquet")
full = pd.read_parquet(tm, columns=["variant_key", "ENSG", "beta_alt"]).merge(
    pd.read_csv(os.path.join(E, "splits/variant_splits.tsv"), sep="\t")[["variant_key", "block", "ld_fold"]],
    on="variant_key", how="inner")
ftr = full[(full["ld_fold"] != j["primary_held_out_fold"]) & full["ENSG"].notna()]
gsum = ftr.groupby("ENSG")["beta_alt"].agg(["sum", "count"])
xx = ftr.join(gsum, on="ENSG")
gbk = ftr.groupby(["ENSG", "block"])["beta_alt"].agg(["sum", "count"]).rename(columns={"sum": "bsum", "count": "bcount"})
xx = xx.join(gbk, on=["ENSG", "block"])
loo = (xx["sum"] - xx["bsum"]) / (xx["count"] - xx["bcount"])
ok = loo.notna() & np.isfinite(loo)
r_in = np.corrcoef(xx["sum"] / xx["count"], xx["beta_alt"])[0, 1]
r_out = np.corrcoef(loo[ok], xx.loc[ok, "beta_alt"])[0, 1]
add("L-15", tm, "training-fold rows (ld_fold != 4, ENSG notna); LD blocks", f"{len(ftr)}; {ftr['block'].nunique()}")
add("L-15", tm, "rows with gene in another block (held-out r rows)", int(ok.sum()))
add("L-15", tm, "recomputed r block left in; r block held out (recorded 0.8669 / 0.0210)", f"{r_in:.5f}; {r_out:.5f}")
rt = p("Analysis/MASLD_Model_Benchmark/executions/routed_multihead_eval-20260829T230000Z/results.json")
rj = json.load(open(rt))
add("L-15", rt, "n_samples; unit_of_inference; donor_key_exists; cohorts",
    f"{rj['n_samples']}; {rj['unit_of_inference']}; {rj['donor_key_exists']}; {'|'.join(rj['bootstrap_distinct_units_per_cohort'])}")

# ---- L-16 / L-17 Fig 7F: reference for per-study deltas, n, Geneformer exposure
inst7f = p("figures/main/fig7_model_framework/fig7F_liver_encoder_against_the_field.pdf")
cand7f = p("figures/candidates/review-r3-fixes-20261015T112148Z/F7a/fig7F_liver_encoder_against_the_field.pdf")
add("L-16", inst7f, "sha256 installed == F7a 7F", sha(inst7f) == sha(cand7f))
txt7f = subprocess.run(["pdftotext", inst7f, "-"], capture_output=True, text=True).stdout
add("L-16", inst7f, "pdftotext contains 'gain over HVG + PCA'", "gain over HVG + PCA" in " ".join(txt7f.split()))
F7 = "figures/candidates/review-r3-fixes-20261015T112148Z/F7a/source_tables/"
gn = p(F7 + "panel_7f_gain.tsv")
gg = pd.read_csv(gn, sep="\t")
bb = gg[gg["block"].str.startswith("b.")]
add("L-16", gn, "block b block_metric", bb["block_metric"].iloc[0])
add("L-16", gn, "block b unit_and_n", bb["unit_and_n"].iloc[0])
psd = p("Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/paired_deltas_by_held_out_study.tsv")
pd_ = pd.read_csv(psd, sep="\t")
add("L-16", psd, "n rows; unique reference_id; unique metric",
    f"{len(pd_)}; {'|'.join(sorted(pd_['reference_id'].unique()))}; {'|'.join(pd_['metric'].unique())}")
eg = p(F7 + "panel_7f_exposure_grid.tsv")
gr = pd.read_csv(eg, sep="\t")
m = gr.merge(pd_[["model_id", "held_out_study", "oriented_delta"]], on=["model_id", "held_out_study"], how="left")
add("L-16", eg, "max |grid per_study_delta - paired_deltas_by_held_out_study oriented_delta|",
    f"{(m['per_study_delta'] - m['oriented_delta']).abs().max():.2e}")
blk = gr[gr["model_id"] == gr["model_id"].iloc[0]]
add("L-16", eg, "one block: sum donors; sum rows; n studies", f"{blk['donors'].sum()}; {blk['rows'].sum()}; {blk['held_out_study'].nunique()}")
oth = gr[gr["held_out_study"] != "GSE189600"]
add("L-16", eg, "other-six-study deltas: n; n positive; min; max",
    f"{len(oth)}; {(oth['per_study_delta'] > 0).sum()}; {oth['per_study_delta'].min():.6f}; {oth['per_study_delta'].max():.6f}")
gf = gr[gr["model_id"] == "geneformer_v2_316m::linear"]
add("L-17", eg, "Geneformer V2 study_exposure_state counts",
    "; ".join(f"{k} {v}" for k, v in gf["study_exposure_state"].value_counts().items()))
add("L-17", eg, "Geneformer V2 encoder_seen studies", "|".join(sorted(gf.loc[gf["study_exposure_state"] == "encoder_seen", "held_out_study"])))
gc_ = gg[gg["block"].str.startswith("c.") & gg["model_id"].str.startswith("geneformer")]
add("L-17", gn, "block c (clean-only) Geneformer gain values", "|".join(str(v) for v in gc_["gain"]))

# ---- L-18 Fig 7 text: no demonstrated advantage over in-cohort ridge
ev = p("Analysis/MASLD_Model_Benchmark/executions/external-validation-expansion-20260901T190000Z/results.json")
vj = json.load(open(ev))
for k in ("GSE213621__primary_3level_no_control", "GSE276114__primary_3level"):
    rr = vj["results"][k]
    add("L-18", ev, f"{k} vs_native_paired_delta; ci95; excludes_zero; verdict",
        f"{rr['vs_native_paired_delta']}; {rr['vs_native_paired_delta_ci95']}; {rr['vs_native_excludes_zero']}; {rr['verdict']}")
rc = vj["reference_cohort_gse268273"]
add("L-18", ev, "GSE268273 native_in_cohort; transfer_spearman (no interval)", f"{rc['native_in_cohort']}; {rc['transfer_spearman']:.4f}")

print("edit_id\tsource_path\tcolumn_or_row\tvalue_reread")
for row in rows:
    print("\t".join(x.replace("\t", " ").replace("\n", " ") for x in row))
