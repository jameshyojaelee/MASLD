"""Checker pass over edits_sc.json: apply corrections and write final_sc.json.

Run: env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python build_final_sc.py
"""
import copy
import json
from pathlib import Path

REPO = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = REPO / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z"
src = json.loads((RUN / "staged/edits_sc.json").read_text())
edits = {e["id"]: copy.deepcopy(e) for e in src}

FIG5F_SRC = "figures/candidates/panel-deep-fixes-20261017T035633Z/fig5/source_tables/5F.tsv"
META = "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_three_cohorts/meta_effects.tsv"
H59 = "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/results_primary_family.tsv"


def fix(eid, reason_add, **kw):
    e = edits[eid]
    for k, v in kw.items():
        e[k] = v
    e["reason"] = e["reason"] + " CHECKER: " + reason_add


fix(
    "sc-02",
    "Neighbors are those of the discovery kNN graph (501_run_hotspot.py create_knn_graph; hotspot/utils.py:48-85, lambda 0.9).",
    new=(
        "Averaging the stored per-cell scores over passing cells keeps Ductular injury donor scores at "
        "Spearman 0.995 with the published scores. Across the 117 programs in the 62-donor equal-variance "
        "fit, stage betas from these passing-cell averages correlate with the published betas at Spearman "
        "0.996, and 2 of 117 change sign. Each stored cell score is 90% an average over the cell's neighbors "
        "in the discovery graph, which include failing cells, so this check is not a rescoring and does not "
        "remove the contribution of failing cells."
    ),
)

fix(
    "sc-03",
    "Named the four gene filters so the 2,000-gene local-correlation step is not read as a depth filter; "
    "depth is computed after strip_confounders and filter_detected (501_run_hotspot.py:262-263, :90-94); "
    "no input h5ad has obs n_counts (checked by h5py this session).",
    new=(
        "For hepatocytes, fibroblasts, macrophages and cholangiocytes, Hotspot used its `danb` model on raw "
        "counts, with each cell's summed counts over the genes kept by the biotype, chromosome, gene-list and "
        "1% detection filters below as the depth term. For T cells, the input matrix held log-normalized "
        "values (log1p of counts per 10,000), and the `danb` model was applied to those values instead of the "
        "integer counts stored in the same file, with the per-cell sum of the log values over the same genes "
        "as the depth term."
    ),
)

fix(
    "sc-04",
    "A corrected integer-scale T-cell score was computed in the diagnostic (complete_raw_vs_corrected_score_r), "
    "so 'no corrected T-cell score' is limited to the native scale; r is Pearson against stored scores "
    "(05_diagnose_tcell_transport.py:116).",
    new=(
        "The 11 prespecified T-cell programs were labeled untestable in this sensitivity. Their native scores "
        "were computed on log-normalized values and can be reconstructed from them (Pearson r with the stored "
        "scores of at least 0.998 for all 11 programs), but the decontX-corrected counts were not transformed "
        "to that scale, so no corrected T-cell score on the native scale was computed."
    ),
)

fix(
    "sc-05",
    "PRESPEC_H59.md:11-21 also lists a two-dataset inverse-variance combination of the focal slopes and says "
    "the H59 answer for the two programs was 'largely predictable'; disclosed. ols_numerator_share is computed "
    "on the 55 fitted donors of the 60-donor set (02_h59_outcome.R:72), not the 62-donor fit, so the "
    "comparison now says 'an equal-variance fit of the same donors'. 'Two named programs' is introduced only "
    "at :357, so the programs are named.",
    new=(
        "The 62-donor equal-variance OLS fit, with BH over 117 programs and an HC3 sensitivity, is a "
        "supplementary sensitivity (Supplementary Note 8). The per-dataset-variance model was specified after "
        "the per-dataset slopes of Ductular injury and Stromal ECM, an inverse-variance combination of those "
        "slopes and the equal-variance results had been seen, and before it was fitted to any program; those "
        "values made its result for these two programs largely predictable. It changes which dataset carries "
        "the estimate: in this model GSE202379 carries 94% (Ductular injury) and 91% (Stromal ECM) of the "
        "information about the stage slope, whereas in an equal-variance fit of the same donors GSE244832 "
        "supplies 84% and 94% of the slope numerator."
    ),
    number_sources=src[4]["number_sources"] + [
        "Analysis/Hypotheses/H59-20261015T112432Z/scripts/02_h59_outcome.R:72 (sxy share over the two fitted datasets)",
    ],
)

fix(
    "sc-06",
    "The GSE213621 control is described in the next paragraph (METHODS.md:405), so a forward pointer is added.",
    new=(
        "No numbered stop rule applied. The program-group permutation control, which failed in GSE213621 "
        "(below), was not run in these two cohorts. A descriptive composition arm"
    ),
)

fix(
    "sc-07",
    "GSA is defined in METHODS only as 'GSA-Human' (key resources table, METHODS.md:24).",
    new=(
        "so 4 of the 32 retained donors are probably controls that cannot be identified. The three excluded "
        "arrays (individuals 27, 28 and 29 in the GSA-Human record) are the only spatial individuals that also "
        "have single-cell data there, and the source study's table lists exactly three spatial donors with "
        "single-cell data, all controls; no source key confirms this join."
    ),
)

fix(
    "sc-08",
    "Rewrote the inference sentence so the arithmetic is visible (1,232,318 printed at fig4.md:3; GSE136103 "
    "61,861 in obs/dataset of scalesc_human_annotated_celltypist.h5ad, read this session). Seed 42 holds for "
    "both builds (03a_export_posterior_only.py:59 random_state=42; 03a_prepare_c2l_reference.py RandomState(42), "
    "file dated 2026-04-07 with the April reference). Training settings for the April build: 250 epochs confirmed "
    "(model.pt history_ elbo_train length 250); 1,000 posterior samples hard-coded in both scripts; batch size and "
    "learning rate come from config/spatial_params.yaml (2500, 0.002), last modified 2026-07-13 [Unverified that "
    "they were unchanged in April].",
    new=(
        "The two Visium sources used different reference builds. For GSE192741, the reference came from an "
        "earlier build of the atlas (1,170,457 cells from six datasets), subsampled to 150,000 cells (seed 42). "
        "Its training log does not name the six datasets. The current atlas has 61,861 more cells (1,232,318), "
        "exactly the number of GSE136103 cells, so GSE136103 was probably the dataset left out. For Vu et al., "
        "the reference came from a later build that includes all seven datasets, with at most 9,375 cells per "
        "cell type (115,562 cells; seed 42). Training used genes detected in at least ten cells, 250 epochs, "
        "batch size 2,500, learning rate 0.002 and 1,000 posterior samples."
    ),
    number_sources=[
        "Analysis/Spatial/logs/c2l_reference_14280337.out:8-14 (1170457 cells, 6 batches, 150000; max_epochs=250)",
        "MASLD_library_design/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad obs/dataset (1,232,318 cells; GSE136103 61,861; read this session); docs/manuscript/draft/fig4.md:3 (1,232,318)",
        "Analysis/Spatial/results/cell2location/reference_model/reference_posterior.h5ad obs (115,562 cells, 7 datasets, max 9,375 per cell_type; read this session)",
        "Analysis/Spatial/results/cell2location/reference_model/model.pt attr_dict history_ elbo_train length 250 (April build)",
        "Analysis/Spatial/scripts/03a_export_posterior_only.py:59 (random_state=42), :75 (min_cells=10), :96 (num_samples 1000); 03a_prepare_c2l_reference.py RandomState(42), :130, :184",
        "Analysis/Spatial/config/spatial_params.yaml:23-25 (max_epochs 250, batch_size 2500, lr 0.002)",
        "Analysis/Spatial/scripts/33_vu_cell2location_deconvolution.py:40",
    ],
)

fix(
    "sc-09",
    "The relaxation widens the expression and detection bins together (04_spatial_projection.py:389-398); "
    "the prior sentence already states hepatocyte matching, so the repeated clause is replaced by the "
    "lineage mismatch.",
    new=(
        "Expression, detection and correlation bins were relaxed by up to two bins when fewer than ten "
        "candidate genes matched; in the four Figure 5F runs this applied to one Vu et al. Stromal ECM gene "
        "(C7), whose expression and detection bins were relaxed by one. Control genes were drawn independently "
        "of one another, so null sets do not reproduce the gene-gene co-expression of a fixed program. Both "
        "programs score highest outside hepatocytes (Stromal ECM in fibroblasts and Ductular injury in "
        "cholangiocytes), and matching did not use correlation with the abundance of those lineages. The null "
        "therefore asks whether a program score is more autocorrelated than a set of independently drawn genes "
        "with matched per-gene properties, not whether it exceeds other co-expressed programs."
    ),
)

fix(
    "sc-10",
    "'pooled null' read as a null distribution; changed to 'pooled estimates'. The OLS numerator share is "
    "for the 55 fitted donors, so 'an equal-variance fit of the same donors'.",
    new=(
        "Each pooled beta is therefore a precision-weighted average, not one stage effect. GSE202379 carries "
        "94% of the information for Ductular injury and 91% for Stromal ECM, so the pooled estimates mainly "
        "reflect GSE202379. In GSE244832 alone (18 donors), the slopes have nominal *p* = 0.011 for Ductular "
        "injury and 0.0089 for Stromal ECM, and an equal-variance fit of the same donors draws most of its "
        "slope from this dataset. We chose the per-dataset-variance model after these slopes and the earlier "
        "equal-variance results were known, and the change moved most of the weight from GSE244832 to "
        "GSE202379 (Methods). Across the 117 programs, the two dataset slopes have opposite signs for 60, and "
        "none has between-dataset heterogeneity at BH q < 0.05 (Cochran Q; smallest q 0.252)."
    ),
    number_sources=src[9]["number_sources"] + [
        H59 + " ols_numerator_share_GSE244832 (0.842632 module 20; 0.940144 module 8)",
    ],
)

fix(
    "sc-12",
    "The 1.58 gap and the 34/12 CI counts are for the arm without composition adjustment (composition arm "
    "min ratio 2.52 at rank 2; CI counts 23/5), so the arm is named. 'Ranks programs largely by agreement, "
    "not by effect size' is false over all 224 rows (Spearman |t| vs |beta| 0.79, vs Q -0.57); replaced by "
    "the verified top-of-ranking statement. Width-Q link verified: corr(log SE, log sqrt(Q/2)) 0.93; "
    "Spearman(CI width, Q) 0.89.",
    new=(
        "No program met BH q < 0.05 on either axis (smallest q 0.079 for fibrosis and 0.159 for NAS), with or "
        "without composition adjustment under either set of cell-type fractions. With three cohorts, the "
        "Knapp-Hartung test uses t on 2 degrees of freedom, so it has low power. Without composition "
        "adjustment, the nearest BH rejection missed by a factor of 1.6 in *p*: four rejections needed a "
        "fourth-smallest *p* of 0.00089 or less, and it was 0.0014. In the same analysis, Knapp-Hartung 95% "
        "intervals exclude zero for 34 fibrosis and 12 NAS programs. These intervals are not adjusted for the "
        "224 tests, and with three cohorts their width depends mainly on each program's between-cohort Q. The "
        "smallest *p* values therefore go to programs whose cohort estimates agree closely, not necessarily to "
        "the largest effects: all 12 of the largest |t| values have Q of 0.49 or less, against an expected "
        "value of 2."
    ),
    number_sources=[
        META + " primary arm: min_k p_(k)/(k*0.05/224) = 1.583 at k = 4; p_(4) = 0.0014138; 4*0.05/224 = 0.000893 (recomputed this session)",
        META + " primary arm: CI excludes 0 for 34/112 fibrosis and 12/112 NAS (recomputed)",
        META + " primary arm: 12 largest |statistic_meta| rows have q_heterogeneity <= 0.490787; expected Q = k - 1 = 2 (recomputed)",
    ],
)

fix(
    "sc-13",
    "'Calls appear only with' dropped 'only': the unprinted common-effect approximation also gives calls "
    "(C10 critic, 45/30 and 24/37).",
    new=(
        "Figure 4E therefore does not resolve whether bulk program scores track linear fibrosis category or "
        "NAS. Calls appear with z-based inference (Figure 4E legend) or with GSE174478 added and the untruncated "
        "Knapp-Hartung variance, and none holds across these analysis choices."
    ),
)

fix(
    "sc-15",
    "Anchor extended to keep the Figure S4 reference on the original sentence; 'over QC-passing cells only "
    "gives' was ambiguous.",
    old="or with a 50-cell lineage gate (q = 0.0765) (Figure S4 [panel letter to be assigned]).",
    new=(
        "or with a 50-cell lineage gate (q = 0.0765) (Figure S4 [panel letter to be assigned]). When the stored "
        "cell scores are averaged over QC-passing cells alone, the call holds (q = 0.0180), but adding a 50- or "
        "20-cell lineage gate to those cells removes it (q = 0.0616 and 0.0889)."
    ),
)

fix(
    "sc-18",
    "'followed cell composition' replaced by the literal regression statement.",
    new=(
        "cell2location abundances and is plotted in Figure 5F, gave the same evidence\n"
        "states with smaller excesses (ECM/IGFBP7 0.040 and 0.050; ductular-injury/BICC1\n"
        "0.015 and 0.003). It removed 71% and 66% of the ECM/IGFBP7 excess, so most of\n"
        "the prespecified excess is accounted for by the estimated abundances of the\n"
        "other cell types. The two sources were deconvolved against different atlas\n"
        "reference builds (Methods)."
    ),
)

fix(
    "sc-19",
    "Moved from the raw-residual sentence (where it split 'The raw statistic therefore' from its premise) to "
    "the donor-level sentence. Exact Ductular injury sign-flip p is 0.1250 in the Fig 5F source table "
    "(panel prints .12 by round-half-even), so 0.125 is printed.",
    old=(
        "At the donor level, ECM/IGFBP7 excess was positive in 3 of 4 donors (exact\n"
        "one-sided sign-flip *p* = 0.19)."
    ),
    new=(
        "At the donor level, excess over the matched null was positive in 3 of 4 donors\n"
        "for ECM/IGFBP7 (exact one-sided sign-flip *p* = 0.19) and for\n"
        "ductular-injury/BICC1 (*p* = 0.125); ductular-injury/BICC1 excess was also\n"
        "positive in 7 of 10 Vu et al. arrays. Figure 5F shows these excess counts, not\n"
        "the raw residual counts above."
    ),
    number_sources=[
        FIG5F_SRC + " pooled rows: Stromal ECM GSE192741 3/4, signflip_p 0.1875; Ductular injury GSE192741 3/4, signflip_p 0.1250; Ductular injury Vu 7/10",
        "figures/main/fig5_molecular_context/panels/fig5f_spatial_program_calibration.pdf (pdftotext: Ductular injury '3/4 donors > 0, P=.12', '7/10 arrays > 0')",
    ],
)

fix(
    "sc-22",
    "Reordered so 'Ductular-injury/BICC1 therefore routes first...' (next sentence) follows its premise.",
    new=(
        "Within the sampled sections, the ECM/IGFBP7 score is more spatially\n"
        "autocorrelated than matched control genes. This is not a donor-level result,\n"
        "and co-expressed programs often exceed the same null. Ductular-injury/BICC1 is\n"
        "spatially indeterminate."
    ),
)

fix(
    "sc-24",
    "Appositive 'matched control genes, a null that' did not parse; restricted 'only' to the two programs.",
    new=(
        "In spatial transcriptomics, of the two programs only ECM/IGFBP7 is more spatially clustered than "
        "matched control genes, a comparison that co-expressed programs often also pass. Ductular-injury/BICC1 "
        "is indeterminate (Figure 5F)."
    ),
)

fix(
    "sc-25",
    "Figure 4E plots fibrosis only; NAS intervals are in Figure S4W.",
    new=(
        "with or without composition adjustment. The intervals are not adjusted for the 224 tests, and with "
        "three cohorts their width depends mainly on between-cohort agreement (Q), so 34 of the 112 fibrosis "
        "intervals (and 12 NAS intervals, Figure S4W) exclude 0 without a BH call. With z-based REML inference "
        "instead, 37 fibrosis and 14 NAS associations reach q < 0.05."
    ),
)

fix(
    "sc-26",
    "Ductular injury sign-flip p printed as the exact 0.125 from the Fig 5F source table.",
    new=(
        "The gene null holds donors and arrays fixed and draws control genes independently, matched on "
        "correlation with hepatocyte abundance, so it is not a donor-level test and is not matched for program "
        "co-expression. Donor and array counts are of units with positive excess over the matched null: "
        "Stromal ECM excess is positive in 3 of 4 Guilliams et al. donors (exact one-sided sign-flip p = 0.19; "
        "the smallest possible p with 4 donors is 0.0625) and Ductular injury excess in 3 of 4 (p = 0.125). "
        "Spot scores use the measured member genes, which carry 55% (Guilliams et al.) and 53% (Vu et al.) of "
        "the Stromal ECM L1 weight and 41% of the Ductular injury weight. IGFBP7 carries 2.47% of the full "
        "Stromal ECM weight and BICC1 5.44% of the full Ductular injury weight, so the result concerns program "
        "scores, not either gene. The two sources used different cell2location reference builds (Methods)."
    ),
    number_sources=[
        "Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-sensitivity-candidate-2026-08-11-r2/full_composition_program_results.tsv retained_l1_weight (0.547515, 0.528185, 0.412663)",
        FIG5F_SRC + " pooled rows signflip_p 0.1875 (Stromal ECM), 0.1250 (Ductular injury), n_units_positive 3/4",
    ],
)

final = [edits[e["id"]] for e in src]

# ---- checks ----
problems = []
by_file = {}
for e in final:
    text = (REPO / e["file"]).read_text()
    n = text.count(e["old"])
    if n != 1:
        problems.append(f"{e['id']}: old occurs {n} times")
    for bad in ("—", "–", "TREAT"):
        if bad in e["new"]:
            problems.append(f"{e['id']}: contains {bad!r}")
    i = text.find(e["old"])
    by_file.setdefault(e["file"], []).append((i, i + len(e["old"]), "sc:" + e["id"]))

# cross-writer overlap with every other staged list present now
for f in sorted((RUN / "staged").glob("*.json")):
    if f.name in ("edits_sc.json", "final_sc.json"):
        continue
    d = json.loads(f.read_text())
    rows = d if isinstance(d, list) else d.get("final_edits", d.get("edits", []))
    for e in rows:
        if e.get("file") in by_file:
            text = (REPO / e["file"]).read_text()
            i = text.find(e["old"])
            if i >= 0:
                by_file[e["file"]].append((i, i + len(e["old"]), f"{f.name}:{e['id']}"))
for rel, spans in by_file.items():
    spans.sort()
    for (a0, a1, ia), (b0, b1, ib) in zip(spans, spans[1:]):
        if b0 < a1 and (ia.startswith("sc:") or ib.startswith("sc:")):
            problems.append(f"overlap in {rel}: {ia} and {ib}")

print("\n".join(problems) if problems else "checks passed")
CHECKER_PROBLEMS = [
    {"id": "sc-12", "kind": "overstated claim", "evidence": "Over all 224 primary-arm rows of meta_effects.tsv, Spearman(|t|, |beta|) is 0.79 and Spearman(|t|, Q) is -0.57, so 'ranks programs largely by agreement, not by effect size' was false; it was replaced by the verified statement that all 12 of the largest |t| have Q <= 0.49 (expected 2). The 1.6 gap and the 34/12 interval counts hold only for the arm without composition adjustment (composition arm: 2.52, 23/5), so that arm is now named."},
    {"id": "sc-04", "kind": "false statement", "evidence": "The original text said 'no corrected T-cell score was computed', but tcell_transport_diagnostic.tsv has complete_raw_vs_corrected_score_r, which is an integer-scale corrected score. The text is now limited to the native scale. The r is Pearson (05_diagnose_tcell_transport.py)."},
    {"id": "sc-05/sc-10", "kind": "wrong comparator", "evidence": "ols_numerator_share_* comes from the 55 fitted donors in the 60-donor set (02_h59_outcome.R:72), not from the 62-donor fit named in the same sentence. The comparator is now 'an equal-variance fit of the same donors'."},
    {"id": "sc-05", "kind": "incomplete disclosure", "evidence": "PRESPEC_H59.md:11-21 says the inverse-variance combination of the two focal slopes was also seen and that the answer was 'largely predictable'. Both are now disclosed."},
    {"id": "sc-19/sc-26", "kind": "rounding", "evidence": "5F.tsv has signflip_p 0.1250 for Ductular injury in GSE192741, and the panel prints P=.12 because of round-half-even. The text now prints 0.125. The panel label in figures/main/fig5 still reads .12; it should read .125 or 0.13."},
    {"id": "sc-19", "kind": "flow", "evidence": "The insertion fell between the 4/4 and 9/10 raw counts and 'The raw statistic therefore does not separate the two programs', so the premise and its conclusion were split. It was moved to the donor-level sentence (fig5.md:102-103)."},
    {"id": "sc-08", "kind": "partially config-sourced", "evidence": "For the April (Vu) build, model.pt and both scripts confirm 250 epochs, min_cells 10, 1,000 samples and seed 42. Batch size 2,500 and lr 0.002 for that build come only from spatial_params.yaml (modified 2026-07-13) [Unverified]. Not checked: whether those values were unchanged in April."},
    {"id": "sc-03", "kind": "unchecked downstream", "evidence": "Not checked: whether the T-cell danb-on-log-values input reaches any printed Fig 5 count (31 and 95). The code default HOTSPOT_MAX_GENES is 5000, while METHODS says 2,000; the T-cell log shows HOTSPOT_MAX_GENES=2000, and the other four logs were not read."},
    {"id": "METHODS gene list", "kind": "possible method misstatement", "evidence": "METHODS says the X-escape gene list applies only to the four non-T inputs, but the current 501_run_hotspot.py strip_confounders applies it to every input. Not checked: whether the code at T-cell run time (job 16236060) matched."},
    {"id": "PAPER.md spatial", "kind": "uncorrected authority", "evidence": "The PAPER.md spatial claim (lines ~314-318, 384-406, 719) still says ECM/IGFBP7 is spatially supported without the co-expression caveat. Left for the AN4 114-program spatial ranking, which is under the hold rule."},
    {"id": "RESULTS.md:905", "kind": "different wording", "evidence": "RESULTS.md says the T-cell substrate 'cannot be reproduced from complete-atlas integer counts'. This does not directly contradict sc-04, which says the scores are reconstructed from log values. Not edited because three other writers edit RESULTS.md."},
    {"id": "sc-15", "kind": "anchor contains placeholder", "evidence": "The new anchor includes 'Figure S4 [panel letter to be assigned]'. A later main_*.json placeholder fill on supplementary_notes.md would overlap. No staged file overlaps it now."},
    {"id": "MANUSCRIPT_r8_assembled.md", "kind": "release copy", "evidence": "The r8 assembled copy (lines 188, 301, 327, 437, 670, 1507, 1521) has the pre-edit text. It is a named release and was not edited."},
]
SUMMARY = ("All 27 sc edits kept, none dropped, none added. 18 rewritten: sc-02, 03, 04, 05, 06, 07, 08, 09, 10, 12, 13, 15, 18, 19, 22, 24, 25, 26. "
           "Substantive fixes: sc-12 claimed the KH test ranks by agreement rather than effect size, which is false over 224 rows; "
           "sc-04 said no corrected T-cell score existed; sc-05/sc-10 compared GLS shares with an OLS fit on different donors; "
           "the Ductular sign-flip p is 0.125, not 0.12; sc-19 was moved to the donor-level sentence. Each old occurs once; "
           "no overlap with edits_bulk, final_bulk, edits_cat or edits_gen; no en or em dashes.")
out = {"final_edits": final, "dropped": [], "added_ids": [], "problems": CHECKER_PROBLEMS, "summary": SUMMARY}
(RUN / "staged/final_sc.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
