"""Stage the lead's text edits (ids T-xx). Every number placed in a `new`
string is read from its producing file in this script. Writes
work/lead_edits.json and work/lead_reread.tsv. Reads only; edits no live file."""
import csv
import json
import os

import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WORK = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-QUALITY-20261017T020832Z/text_fixes-20261017T025625Z/work")
os.chdir(ROOT)

R9D = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet"
H55 = "Analysis/Hypotheses/H55-20261015T112558Z/run-20261015T114811Z/outputs"
H59 = "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/results_primary_family.tsv"
H81 = "Analysis/Hypotheses/H81-20261015T125233Z/run-20261015T141242Z/summary/O1_detection_limit_by_decile.tsv"
FIG1C = "figures/main/fig1_atlas_overview/fig1c_trait_scope_source.tsv"
FIG3G = "figures/candidates/fig3g-working-layout-774-20261012T011006Z/source_tables/fig3g_raw_treat_interface.tsv"
CLCE = "figures/main/fig4_singlecell_programs/source_tables/fig34_scientific_corrections_20260930/content_lineage_cohort_effects.tsv"
CFT = "RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T141625Z/PRIMARY_transfer_verdict.tsv"
PXD = "data/PXD051911/meta_data.txt"
SIG774 = "RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261011T152106Z/bulk/run_summary.tsv"
STUDY = "Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/study_metrics.tsv"
P7C = "Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/tables/panel_7c_exposure.tsv"
ROSTER4H = "figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/roster_census.tsv"
SUM4H = "figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/summary_stage_vs_continuum_by_roster.tsv"
INSTALL = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave4/INSTALL_PLAN.tsv"
OR774 = "figures/candidates/fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_matched_treat_summary.tsv"
OR929 = "figures/candidates/fig34-scientific-corrections-20260929-v1/figures_v2/fig3g_matched_uncertainty.tsv"
RED774 = "RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/hac-system-stress-corrected774-20261009T171017Z/tables/system_redundancy.tsv"
RED844 = "RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/hac-system-stress-20260827T125516Z/tables/system_redundancy.tsv"
RESULTS_MD = "docs/RESULTS.md"

reread_rows = []


def note(eid, path, what, value):
    reread_rows.append([eid, path, what, str(value)])


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------- re-reads ----------
r9 = pd.read_parquet(R9D, columns=[
    "symbol", "primary_evidence_class", "next_experiment_rule_id",
    "r8_genetics_call_state", "gen_coloc_best_susie_pp4", "gen_coloc_best_susie_study",
    "gen_coloc_best_abf_pp4", "gen_coloc_best_abf_study", "coding_cs_max_pip",
    "coding_cs_state", "caqtl_region_state", "gtex_postmortem_fat_state"])
gnmt = r9[r9.symbol == "GNMT"].iloc[0]
gnmt_abf = f"{gnmt.gen_coloc_best_abf_pp4:.3f}"
gnmt_susie = f"{gnmt.gen_coloc_best_susie_pp4:.3f}"
for k in ["gen_coloc_best_abf_pp4", "gen_coloc_best_abf_study", "gen_coloc_best_susie_pp4",
          "gen_coloc_best_susie_study", "r8_genetics_call_state", "primary_evidence_class",
          "next_experiment_rule_id"]:
    for eid in ["T-02", "T-03", "T-04"]:
        note(eid, R9D, f"GNMT {k}", gnmt[k])
assert gnmt.r8_genetics_call_state == "indeterminate"
assert gnmt.gen_coloc_best_abf_study == "2021_34841290_NAFLD_EUR"

f1c = pd.read_csv(FIG1C, sep="\t")
n_direct = int(f1c.loc[f1c.genetic_trait_scope.isin(["direct_disease", "both"]), "n_genes"].sum())
note("T-01", FIG1C, "n_genes direct_disease + both", n_direct)

h59 = pd.read_csv(H59, sep="\t")
h59_min_q = f"{h59.perm_q.min():.3f}"
for eid in ["T-05", "T-06"]:
    note(eid, H59, "n programs; sum(selected_H59); min(perm_q)",
         f"{len(h59)}; {int(h59.selected_H59.sum())}; {h59.perm_q.min()}")
assert len(h59) == 117 and int(h59.selected_H59.sum()) == 0

fam = pd.read_csv(os.path.join(H55, "family_slopes.tsv"), sep="\t")
rho = fam[["bulk_logFC", "slope"]].corr(method="spearman").iloc[0, 1]
col = pd.read_csv(os.path.join(H55, "catalog_column_gtex_postmortem_fat.tsv"), sep="\t")
measurable = col.gtex_postmortem_fat_q.notna().sum()
n_q05 = int((col.gtex_postmortem_fat_q < 0.05).sum())
note("T-07", H55 + "/family_slopes.tsv", "n genes; Spearman(bulk_logFC, slope)", f"{len(fam)}; {rho:.6f}")
note("T-07", H55 + "/catalog_column_gtex_postmortem_fat.tsv", "genes with q; q<0.05", f"{measurable}; {n_q05}")
r9_meas = int(r9.gtex_postmortem_fat_state.isin(["not_detected", "associated_up_agrees", "associated_down_agrees", "associated_down_disagrees"]).sum())
r9_assoc = int(r9.gtex_postmortem_fat_state.str.startswith("associated").sum())
note("T-07", R9D, "gtex_postmortem_fat_state measurable; associated", f"{r9_meas}; {r9_assoc}")
assert len(fam) == measurable == r9_meas and n_q05 == r9_assoc

coding31 = r9[r9.coding_cs_state == "fine_mapped_coding_pip_ge_0.5"]
note("T-08", R9D, "coding_cs_state fine_mapped_coding_pip_ge_0.5 n; classes",
     f"{len(coding31)}; {coding31.primary_evidence_class.value_counts().to_dict()}")
susie436 = r9[r9.gen_coloc_best_susie_pp4 > 0.5]
note("T-09", R9D, "SuSiE>0.5 genes; caqtl 'colocalized accessible region'",
     f"{len(susie436)}; {(susie436.caqtl_region_state == 'colocalized accessible region').sum()}")
for eid in ["T-08", "T-09", "T-10", "T-11", "T-12", "T-13"]:
    note(eid, RESULTS_MD, "adoption line", "Catalog r9d was adopted on 2026-10-16" in open(RESULTS_MD).read())
pn = r9[r9.symbol == "PNPLA3"].iloc[0].coding_cs_max_pip
mt = r9[r9.symbol == "MTARC1"].iloc[0].coding_cs_max_pip
note("T-12", R9D, "coding_cs_max_pip PNPLA3; MTARC1", f"{pn:.6f}; {mt:.6f}")
note("T-13", R9D, "column gtex_postmortem_fat_state present", "gtex_postmortem_fat_state" in r9.columns)

h81 = pd.read_csv(H81, sep="\t", dtype={"detection_limit_log2FC": str})
num = h81[h81.detection_limit_log2FC != "above 1.5"].copy()
num["v"] = num.detection_limit_log2FC.astype(float)
best = num.loc[num.v.idxmin()]
dl_min, dl_min_dec = best.detection_limit_log2FC, int(best.decile)
dl_top = h81.loc[h81.decile == 10, "detection_limit_log2FC"].iloc[0]
dl_bottom = h81.loc[h81.decile == 1, "detection_limit_log2FC"].iloc[0]
for eid in ["T-16", "T-17"]:
    note(eid, H81, "min limit (decile); decile 10; decile 1", f"{dl_min} ({dl_min_dec}); {dl_top}; {dl_bottom}")
assert dl_bottom == "above 1.5"

g3 = pd.read_csv(FIG3G, sep="\t")
rho3g = g3.logFC.abs().corr(g3.susie_pp4, method="spearman")
note("T-20", FIG3G, "n pairs; Spearman(|logFC|, susie_pp4)", f"{len(g3)}; {rho3g:.6f}")
rho3g_txt = f"{rho3g:.4f}"

cl = pd.read_csv(CLCE, sep="\t")
sel = cl[(cl.dataset == "GSE202379") & (cl.minimum_cells == 50) & (cl.annotation_filter == "all")
         & (cl.arm == "ambient_corrected_original_exposure") & (cl.contrast == "secondary_steatohepatitis_vs_healthy")]
nref = {r.program: (int(r.n_reference), r.lineage, r.failure_reason) for r in sel.itertuples()}
note("T-21", CLCE, "GSE202379 50-cell healthy-reference donors", nref)
n_chol, n_fib = nref["Ductular injury"][0], nref["Stromal ECM"][0]
oth = cl[(cl.dataset == "GSE244832") & (cl.minimum_cells == 50) & (cl.annotation_filter == "all")
         & (cl.arm == "ambient_corrected_original_exposure") & (cl.contrast == "secondary_steatohepatitis_vs_healthy")]
note("T-21", CLCE, "GSE244832 50-cell healthy-reference donors", {r.program: int(r.n_reference) for r in oth.itertuples()})

cft = pd.read_csv(CFT, sep="\t")
note("T-23", CFT, "cohort; observed_rho; p_holm_2test",
     "; ".join(f"{r.cohort} {r.observed_rho:.3f} {r.p_holm_2test:.3f}" for r in cft.itertuples()))
assert (cft.p_holm_2test > 0.05).all()

meta = pd.read_csv(PXD, sep="\t", keep_default_na=False)
liv = meta[~meta.liver_proteomics_filename.isin(["NA", ""])]
n_liver = liv.patient_name.nunique()
n_ctrl = int((liv.saf_diagnosis == "No_MASLD").sum())
n_masld = int(liv.saf_diagnosis.isin(["MASL", "MASH"]).sum())
note("T-24", PXD, "liver participants; MASL+MASH; No_MASLD", f"{n_liver}; {n_masld}; {n_ctrl}")
assert n_liver == n_masld + n_ctrl

sig = pd.read_csv(SIG774, sep="\t").set_index("metric").value
n_pub, n_obs = int(sig["signature_published"]), int(sig["signature_observed"])
note("T-26", SIG774, "signature_published; signature_observed", f"{n_pub}; {n_obs}")

sm = pd.read_csv(STUDY, sep="\t")
tf = sm[sm.model_id == "transcriptformer_tf_sapiens::two_layer_mlp"]
n_cells, n_donors, n_studies = int(tf.rows.sum()), int(tf.donors.sum()), tf.held_out_study.nunique()
note("T-27", STUDY, "TF-Sapiens two_layer_mlp: sum rows; sum donors; studies", f"{n_cells}; {n_donors}; {n_studies}")
p7c = pd.read_csv(P7C, sep="\t", comment=None)
p7c = p7c[p7c.block != "block"]
tf7c = p7c[p7c.model.str.contains("transcriptformer", case=False)]
note("T-27", P7C, "7C TranscriptFormer unit_and_n; score", f"{tf7c.unit_and_n.iloc[0]}; {float(tf7c.score.iloc[0]):.4f}")
words = {7: "seven"}

ros = pd.read_csv(ROSTER4H, sep="\t")
dis = ros[ros.group_binary == "Disease"].groupby("dataset").n.sum().to_dict()
n_dis = sum(dis.values())
n_ctl = int(ros[ros.group_binary == "Control"].n.sum())
note("T-28", ROSTER4H, "disease units by dataset; controls", f"{dis}; total {n_dis}; controls {n_ctl}")
inst = pd.read_csv(INSTALL, sep="\t")
g_rows = inst[inst["item"] == "G"]
note("T-28", INSTALL, "item G targets", "; ".join(os.path.basename(t) for t in g_rows.target))
inst_dir = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave4/install-20261016T104914Z"
note("T-28", inst_dir, "install directory exists (2026-10-16)", os.path.isdir(inst_dir))
for src, tgt in zip(g_rows.source, g_rows.target):
    note("T-28", tgt, "sha256 installed == item G source", sha256(src) == sha256(tgt))
s4h = pd.read_csv(SUM4H, sep="\t")
row = s4h[(s4h.level == "systems_43") & (s4h.roster == "disease_only") & (s4h.x == "stage") & (s4h.y == "cont_fixed_projection")].iloc[0]
note("T-29", SUM4H, "disease_only systems_43 stage vs fixed projection rho; same_direction",
     f"{row.spearman_rho:.4f}; {row.same_direction}/{row.n}")

o774 = pd.read_csv(OR774, sep="\t").iloc[0]
o929 = pd.read_csv(OR929, sep="\t").iloc[0]
note("T-30", OR774, "n pairs; paired_or; paired_p", f"{o774.n_expression_matched_pairs}; {o774.paired_or:.5f}; {o774.paired_p:.6f}")
note("T-30", OR929, "n pairs; OR; CI; p", f"{o929.n_pairs}; {o929.paired_OR}; {o929.OR_low:.3f}-{o929.OR_high:.3f}; {o929.nominal_exact_McNemar_p:.3f}")

red7 = pd.read_csv(RED774, sep="\t").set_index("scope")
red8 = pd.read_csv(RED844, sep="\t").set_index("scope")
for nm, r in [("774", red7), ("844", red8)]:
    note("T-31", RED774 if nm == "774" else RED844, f"{nm}: pooled n; eff; pc1; pc1-3; GSE162694 eff; GSE213621 eff",
         f"{r.loc['pooled_centred','n_participants']}; {r.loc['pooled_centred','effective_independent_systems']}; "
         f"{r.loc['pooled_centred','pc1_variance_explained']:.4f}; {r.loc['pooled_centred','pc1_to_pc3_variance_explained']:.4f}; "
         f"{r.loc['GSE162694','effective_independent_systems']}; {r.loc['GSE213621','effective_independent_systems']}")


def pct(x):
    return f"{100 * x:.1f}%"


p7 = red7.loc["pooled_centred"]
p8 = red8.loc["pooled_centred"]

# ---------- edits ----------
E = []


def add(eid, file, old, new, reason, source):
    E.append(dict(id=eid, file=file, old=old, new=new, reason=reason, source=source))


add("T-01", "docs/manuscript/draft/01_intro.md",
    "22 genes have\nany direct-trait colocalization.",
    f"{n_direct} genes have\na direct-trait multi-signal SuSiE colocalization (PP.H4 above 0.5).",
    "Count is multi-signal SuSiE only (ABF gives a larger direct set); match fig2.md wording 'direct-trait SuSiE PP.H4 above 0.5'.",
    f"{FIG1C} (n_genes, rows direct_disease + both)")

add("T-02", "docs/manuscript/draft/01_intro.md",
    "GNMT begins\nas a provisional direct-disease shared-signal locus with an unresolved local\n"
    "regulatory haplotype, then joins MAT1A and CYP2C19 in lower\n"
    "disease-state RNA, cross-sectional fibrosis-stage abundance, and liver protein.\n"
    "This locus-to-state bridge connects the maps without claiming that the\n"
    "three-gene metabolic state is one inherited mechanism.",
    f"GNMT has a\nsingle-signal direct-NAFLD colocalization (ABF PP.H4 {gnmt_abf}) but no multi-signal\n"
    f"SuSiE call (best PP.H4 {gnmt_susie}), so its Catalog genetics call is indeterminate,\n"
    "and its local regulatory haplotype is unresolved. With MAT1A and CYP2C19 it has lower\n"
    "disease-state RNA, cross-sectional fibrosis-stage abundance, and liver protein.\n"
    "The three-gene metabolic state is not claimed to be one inherited mechanism.",
    "GNMT has no multi-signal call and an indeterminate Catalog genetics call; 'provisional direct-disease locus' and 'bridge connects the maps' state a genetic link the Resource does not call.",
    f"{R9D} (GNMT: gen_coloc_best_abf_pp4, gen_coloc_best_abf_study, gen_coloc_best_susie_pp4, r8_genetics_call_state)")

add("T-03", "docs/manuscript/draft/fig5_discussion.md",
    "direct-NAFLD ABF result (PP.H4 0.706), so GNMT also routes to allele-aware resolution of its local regulatory haplotype; this\n"
    "does not extend to MAT1A or CYP2C19.",
    f"single-signal direct-NAFLD ABF result (PP.H4 {gnmt_abf}); this\n"
    f"does not extend to MAT1A or CYP2C19. GNMT has no multi-signal SuSiE call (best PP.H4 {gnmt_susie}), its Catalog genetics call is indeterminate, and its Catalog route is the same perturbation and rescue route; its local regulatory haplotype is unresolved.",
    "GNMT's Catalog route is EXP_STRESS_PERTURB_V1 (knockdown and rescue), not allele-aware resolution; its genetics call is indeterminate.",
    f"{R9D} (GNMT: gen_coloc_best_abf_pp4, gen_coloc_best_susie_pp4, r8_genetics_call_state, next_experiment_rule_id)")

add("T-04", "docs/manuscript/draft/fig5_discussion.md",
    "GNMT is a provisional direct-disease locus that intersects a coordinated\n"
    "GNMT–MAT1A–CYP2C19 loss across pooled RNA, cross-sectional fibrosis contrasts,\n"
    "and adjusted liver protein.",
    f"GNMT has a single-signal direct-NAFLD colocalization (ABF PP.H4 {gnmt_abf}) but no\n"
    "multi-signal call, and with MAT1A and CYP2C19 it is lower in pooled RNA, cross-sectional fibrosis contrasts,\n"
    "and adjusted liver protein.",
    "'Direct-disease locus' rests on one single-signal ABF result; the multi-signal call is indeterminate. 'Loss' replaced by 'lower' (cross-sectional).",
    f"{R9D} (GNMT: gen_coloc_best_abf_pp4, gen_coloc_best_susie_pp4, r8_genetics_call_state)")

add("T-05", "docs/manuscript/draft/abstract.md",
    "none passed the 117-program family test for a cross-sectional stage\nassociation.",
    f"none passed within-dataset permutation BH over the {len(h59)} programs for a\ncross-sectional stage association (smallest q {h59_min_q}).",
    "Name the correction and family (RESULTS: within-dataset permutation BH over all 117 programs; smallest q=0.323375).",
    f"{H59} (rows; selected_H59; perm_q)")

add("T-06", "docs/manuscript/draft/01_intro.md",
    "Neither focal program passes the 117-program rule for a\ncross-sectional stage association once dataset-specific variances are\nmodeled.",
    f"Neither focal program passes within-dataset permutation BH over all {len(h59)} programs for a\n"
    f"cross-sectional stage association once dataset-specific variances are\nmodeled (smallest q of any program {h59_min_q}).",
    "Name the correction and family instead of '117-program rule'.",
    f"{H59} (rows; selected_H59; perm_q)")

add("T-07", "docs/manuscript/draft/abstract.md",
    "the 1,002 disease-state genes measurable in GTEx changed with image-measured\n"
    "histologic fat in the same direction as in the biopsy disease-versus-control\n"
    "contrast (Spearman rho 0.638)",
    f"across the {measurable:,} disease-state genes measurable in GTEx, each gene's\n"
    "association with image-measured histologic fat correlated with its biopsy\n"
    f"disease-versus-control log2 fold change (gene-level Spearman rho {rho:.3f}; {n_q05}\n"
    f"of the {measurable:,} genes have a fat slope at BH q < 0.05 over those genes)",
    "The 0.638 is a gene-level correlation; only 326 of 1,002 genes reach BH q<0.05 individually, so the text must not imply all 1,002 changed.",
    f"{H55}/family_slopes.tsv (bulk_logFC, slope: Spearman); {H55}/catalog_column_gtex_postmortem_fat.tsv (gtex_postmortem_fat_q); cross-checked in {R9D} (gtex_postmortem_fat_state)")

add("T-08", "docs/manuscript/draft/fig2.md",
    "A candidate Catalog field,\nnot part of r8, records such members:",
    "A field of Catalog r9d (adopted 2026-10-16),\nnot part of r8, records such members:",
    "H87 coding field is an adopted r9d field (RESULTS 'Catalog r9d fields, adopted 2026-10-16'); counts unchanged.",
    f"{RESULTS_MD} ('Catalog r9d was adopted on 2026-10-16'); counts re-read from {R9D} (coding_cs_state, primary_evidence_class)")

add("T-09a", RESULTS_MD,
    "Two candidate annotation fields, not part of r8, add\nto these records",
    "Two annotation fields of Catalog r9d (adopted 2026-10-16), not part of r8, add\nto these records",
    "Release state: these fields are adopted in r9d (r9d = r9c corrections of the r9 fields plus H79/H81).",
    f"{RESULTS_MD} r9d adoption section; {R9D}")
add("T-09b", RESULTS_MD,
    "A separate\ncandidate field records fine-mapped protein-altering",
    "A separate\nr9d field records fine-mapped protein-altering",
    "Release state: H87 field is adopted in r9d; 31 genes re-read.",
    f"{R9D} (coding_cs_state fine_mapped_coding_pip_ge_0.5)")

add("T-10", "docs/manuscript/draft/fig5_discussion.md",
    "The r8 Gene Catalog is selected for the paper;",
    "Gene Catalog r9d (adopted 2026-10-16), which keeps the r8 classes, is selected for the paper;",
    "Adopted Catalog is r9d; r8 classes unchanged.",
    f"{RESULTS_MD} ('Catalog r9d was adopted on 2026-10-16 ... keeps ... all r8 class memberships')")
add("T-11", "docs/manuscript/draft/fig5_discussion.md",
    "The r9 candidate Catalog adds fields",
    "Catalog r9d (adopted 2026-10-16) adds fields",
    "r9 candidate fields are now adopted r9d fields.",
    f"{RESULTS_MD} r9d adoption section")
add("T-12", "docs/manuscript/draft/fig5_discussion.md",
    "In the r9 candidate coding field each carries a",
    "In the r9d coding field each carries a",
    "Coding field is adopted in r9d; PIP values in the following text re-read and unchanged.",
    f"{R9D} (coding_cs_max_pip PNPLA3, MTARC1)")
add("T-13", RESULTS_MD,
    "`gtex_postmortem_fat` annotation is for the next Catalog build, not r8.",
    "`gtex_postmortem_fat` annotation is a field of Catalog r9d (adopted 2026-10-16), not r8.",
    "Extra (same class as the r9d label fixes): the H55 column is in the adopted r9d table.",
    f"{R9D} (column gtex_postmortem_fat_state)")

add("T-14", RESULTS_MD,
    "The Catalog still shows 0 strict `tested_negative` calls; adopting the\nrule is a user decision",
    "The Catalog `tested_negative` state remains not assessed (no negative-call\nrule is defined); adopting the rule is a user decision",
    "Contradicts the Catalog row 'not assessed' (no negative-call rule; a count of zero would be fixed by construction).",
    f"{RESULTS_MD} Catalog row 'Catalog `tested_negative` calls | not assessed'")
add("T-15", RESULTS_MD,
    "Not adopted; the Catalog count stays 0 |",
    "Not adopted; the Catalog `tested_negative` state stays not assessed |",
    "Same contradiction as T-14, in the observability table.",
    f"{RESULTS_MD} Catalog row 'Catalog `tested_negative` calls | not assessed'")

old_dl = ("the smallest change detected with 80% recall ranged\n"
          "from 0.653 log2FC in the most highly expressed decile to above 1.5 in the\nleast.")
new_dl = (f"the smallest change detected with 80% recall was\n"
          f"{dl_min} log2FC at its lowest (decile {dl_min_dec}), {dl_top} in the most highly expressed decile and\n"
          f"{dl_bottom} in the least.")
add("T-16", "docs/manuscript/draft/fig3.md", old_dl, new_dl,
    "'ranged from 0.653' reads as the lower end; the smallest limit is 0.649 in decile 8.", f"{H81} (decile, detection_limit_log2FC)")
add("T-17", RESULTS_MD, old_dl, new_dl,
    "Same phrasing in the RESULTS H81 paragraph.", f"{H81} (decile, detection_limit_log2FC)")

add("T-18", "docs/manuscript/draft/fig3.md",
    "All three decrease in the pooled disease-state comparison",
    "All three are lower in the pooled disease-state comparison",
    "Cross-sectional comparisons; 'decrease' implies change over time.", "wording only")
add("T-19", "docs/manuscript/draft/fig3.md",
    "which bounds the FDR over all four 24,115-gene families at 5%;",
    "which bounds the FDR over all discoveries in the four 24,115-gene families at 5% but not specifically within the SuSiE subset;",
    "Add the RESULTS qualifier so the 15 genes are not read as 5% FDR-controlled within the SuSiE subset.",
    f"{RESULTS_MD} Figure 3F paragraph ('bounds the FDR over all discoveries in the four families at 5% but not specifically within the SuSiE subset')")

add("T-20", "docs/manuscript/05_figure_legends.md",
    "(Spearman rho=0.001, nominal P=0.93)",
    f"(Spearman rho={rho3g_txt}, nominal P=0.93)",
    "Match RESULTS and fig3.md precision (0.0014).", f"{FIG3G} (Spearman of |logFC| with susie_pp4, 3,776 rows)")

add("T-21", "docs/manuscript/draft/fig4.md",
    "because the smaller cohort lacks three\nhealthy donors in these lineages.",
    f"because GSE202379 has fewer than three\nhealthy donors in these lineages ({n_chol} cholangiocyte and {n_fib} fibroblast reference donors at the 50-cell threshold).",
    "The limiting cohort is GSE202379 (the larger single-cell cohort), not 'the smaller cohort'; GSE244832 has 5 healthy reference donors.",
    f"{CLCE} (dataset GSE202379, minimum_cells 50, arm ambient_corrected_original_exposure, contrast secondary_steatohepatitis_vs_healthy: n_reference, failure_reason)")

add("T-22", "docs/manuscript/draft/fig4.md",
    "these are not exact F0–F3\nstages.",
    "these are not exact F0–F4\nstages.",
    "The GSE213621 bins span F0 to F4 (F3–F4 is one bin).", "wording; bins as quoted in the same sentence")

add("T-23", "docs/manuscript/draft/fig4.md",
    "retain matched-background Holm P=0.117 and 0.064.",
    "have matched-background Holm P={:.3f} and {:.3f}, neither below 0.05.".format(
        cft.set_index("cohort").loc["GSE174478", "p_holm_2test"], cft.set_index("cohort").loc["GSE240729", "p_holm_2test"]),
    "'retain' reads as if the checks pass; neither Holm P is below 0.05. Values are from the 2026-08-28 pre-774 run (see needs_user.md).",
    f"{CFT} (cohort, observed_rho, p_holm_2test)")

add("T-24", "docs/manuscript/draft/fig5.md",
    "liver proteome (58 patients).",
    f"liver proteome ({n_liver} participants: {n_masld} MASLD, {n_ctrl} controls).",
    "12 of the 58 are controls; 'patients' is wrong for them.",
    f"{PXD} (rows with liver_proteomics_filename; saf_diagnosis MASL/MASH vs No_MASLD)")

add("T-25", "docs/manuscript/draft/fig5_discussion.md",
    "reported a moderate RNA-protein zonation correlation against single-nucleus RNA",
    "reported an RNA-protein zonation comparison against single-nucleus RNA",
    "'moderate' is not in RESULTS; RESULTS describes the source as reporting a comparison.",
    f"{RESULTS_MD} H63 paragraph ('Weiss et al. reported an RNA-protein zonation comparison against single-nucleus RNA')")

add("T-26", "docs/manuscript/05_figure_legends.md",
    "after removing the 145 signature genes.",
    f"after removing the {n_obs} observed genes of the published {n_pub}-gene signature.",
    "S3S-U are on the 774-unit basis, where 143 of the 145 signature genes are observed (the check's 139 is the 844-record value).",
    f"{SIG774} (signature_published, signature_observed); RESULTS 774-unit rerun paragraph states the same 143")

add("T-27", "docs/manuscript/05_figure_legends.md",
    "its pooled broad cell-class macro F1 of 0.9389\nincludes both",
    f"its pooled broad cell-class macro F1 of {float(tf7c.score.iloc[0]):.4f}\n({n_cells:,} donor-balanced cells from {n_donors} donors in {words[n_studies]} held-out studies) includes both",
    "Add the 7C n; the 7C cell-state rows use the same 50,000-cell evaluation as 7F.",
    f"{STUDY} (transcriptformer_tf_sapiens::two_layer_mlp: rows, donors, held_out_study); {P7C} (unit_and_n, score)")

add("T-28a", RESULTS_MD,
    "current Figure 4G/H installed on 2026-10-14)",
    "Figure 4G/H installed on 2026-10-14 and replaced on 2026-10-16 by the renders on the 367 disease units)",
    "The installed 4G/4H are the disease-only renders (round-4 install item G, 2026-10-16; sha256 match).",
    f"{INSTALL} (item G); sha256 of installed 4G/4H vs item G sources; {ROSTER4H}")
add("T-28b", RESULTS_MD,
    "so the installed Figure 4G and 4H models include those controls as their own stage category and the 4H stage axis reads the control label (review round 3, R05-02; candidate 4G and 4H on the 367 disease units in",
    f"so the earlier Figure 4G and 4H models included those controls as their own stage category and the 4H stage axis read the control label; the installed Figure 4G and 4H (2026-10-16) are the refits on the {n_dis} disease units, {dis['GSE162694']} GSE162694 and {dis['GSE213621']} GSE213621, without the {n_ctl} controls (review round 3, R05-02;",
    "Sentence describes the replaced panels as installed.",
    f"{ROSTER4H} (n by dataset and group); {INSTALL} (item G); sha256 match of installed PDFs")
add("T-29", RESULTS_MD,
    "(printed as 0.73) and 37/43 in the same direction; with the 68 GSE213621\n"
    "controls removed from both axes (367 disease units) the values are 0.631 and\n"
    "32/43 (candidate `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`)",
    "(printed as 0.73 on the earlier panel) and 37/43 in the same direction; with the 68 GSE213621\n"
    f"controls removed from both axes (367 disease units) the values are {row.spearman_rho:.3f} and\n"
    f"{row.same_direction}/{row.n} (installed Figure 4H since 2026-10-16; `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`)",
    "The 367-unit values are the installed panel, not a candidate.",
    f"{SUM4H} (systems_43, disease_only, stage vs cont_fixed_projection); {INSTALL} (item G)")

add("T-30", RESULTS_MD,
    "Figure 3G matched interval-null sensitivity: OR=1.32, exact 95% CI 0.761–2.315,\n  nominal McNemar P=0.358.",
    f"Figure 3G matched interval-null sensitivity, superseded value from this 2026-09-29 candidate ({o929.n_pairs} pairs): OR=1.32, exact 95% CI 0.761–2.315,\n"
    f"  nominal McNemar P=0.358. Current 774-unit value ({o774.n_expression_matched_pairs} pairs): paired OR={o774.paired_or:.3f}, nominal exact McNemar P={o774.paired_p:.3f}\n"
    "  (`figures/candidates/fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_matched_treat_summary.tsv`).",
    "Unlabelled historical value; current value is in RESULTS S3G paragraph.",
    f"{OR929} (n_pairs, paired_OR, OR_low, OR_high, nominal_exact_McNemar_p); {OR774} (n_expression_matched_pairs, paired_or, paired_p)")

add("T-31", RESULTS_MD,
    "The 43 systems are also not 43 independent quantities. Across the pooled\n"
    "evaluation participants their score correlation matrix gives a Li and Ji\n"
    "effective number of 15 independent systems, with 36.5% of variance in the first\n"
    "component and 72.4% in the first three; the per-cohort effective numbers are 12\n"
    "and 14.",
    "The 43 systems are also not 43 independent quantities. Across the "
    f"{int(p7.n_participants)} pooled\n"
    "evaluation units of the 774-unit substrate their score correlation matrix gives a Li and Ji\n"
    f"effective number of {int(p7.effective_independent_systems)} independent systems, with {pct(p7.pc1_variance_explained)} of variance in the first\n"
    f"component and {pct(p7.pc1_to_pc3_variance_explained)} in the first three; the per-cohort effective numbers are "
    f"{int(red7.loc['GSE162694','effective_independent_systems'])}\n"
    f"and {int(red7.loc['GSE213621','effective_independent_systems'])} (`{RED774}`). Superseded 2026-10-09: on the 844-record basis "
    f"({int(p8.n_participants)} units) the values were {int(p8.effective_independent_systems)} systems, "
    f"{pct(p8.pc1_variance_explained)} and {pct(p8.pc1_to_pc3_variance_explained)}.",
    "Paragraph printed the 844-record values without a basis label; replace with the 774-unit values and label the old ones superseded.",
    f"{RED774} and {RED844} (pooled_centred and per-cohort rows)")

with open(os.path.join(WORK, "lead_edits.json"), "w") as f:
    json.dump(E, f, indent=1, ensure_ascii=False)
with open(os.path.join(WORK, "lead_reread.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["edit_id", "source_path", "column_or_row", "value_reread"])
    w.writerows(reread_rows)
print(len(E), "lead edits;", len(reread_rows), "reread rows")
