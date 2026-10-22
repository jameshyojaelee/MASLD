"""Assemble the F3ef repair entries (structured output payload) from text_edits_repair.json."""
import json
B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
RD = B + "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
F1 = B + "figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/"
FR = F1 + "repair/"
TE = json.load(open(RD + "repair/F3ef/text_edits_repair.json"))
first = {e["fid"]: e for e in json.load(open(RD + "repair_inputs/F3ef.json"))["entries"]}
REFIT = "sbatch 7550 (slurm 22256140, COMPLETED 00:04:02, 8 CPU / 64G) running repair/F3ef/refit-20261015T171823Z/refit_cohort_omission.R (31 fits; model copied from scripts/figures/build_pi_stage_extensions.R:142-241); rnaseq python repair/F3ef/scripts/summarize_refit.py"
H82 = "H82 (gate G2-02) wrote RESULT_H82.md at 2026-10-15 17:22Z, after plan_results.json listed it as running. Verdict computed in code: UNTESTABLE (S5, planted-covariate control PC-3 failed); 'No H82 number enters the paper. The review numbers stay review-only.'"

E = []
E.append(dict(fid="R03-03", status="fixed", changed_in_repair=True,
  reproduced=first["R03-03"]["reproduced"] + " Repair: job 7550 roster assert gives 387 units in the three cohorts.",
  what_changes="No number changes. The Fig 3 legend says all panels use the 774-unit substrate; the legend now says all panels draw on it and the 3E/3F stage contrasts use 387 units from three cohorts. Legend, draft and METHODS name the three cohorts and n=387 and say GSE213621 records only binned fibrosis. Candidate 3A (round 1) shows GSE213621 fibrosis as 'Binned' instead of 'Available'.",
  new_value=first["R03-03"]["new_value"], old_value=first["R03-03"]["old_value"],
  panels=first["R03-03"]["panels"], source_tables=first["R03-03"]["source_tables"],
  commands=first["R03-03"]["commands"] + "; repair: old_text re-verified once each in current files (repair/F3ef/scripts/build_text_edits.py)",
  text_edits=TE["R03-03"],
  note="Repair changes: 'exact Kleiner stage' replaced by 'an exact F0–F4 fibrosis stage' (no repository document confirms the scoring system); legend 113 now reads 'All panels draw on the corrected 774-unit substrate; the stage contrasts in E and F use the 387 units ...' (checker wording, avoids implying D uses all 774). 3A belongs to the opening-panel set; the legend 116-117 sentence is true with either the installed or the candidate 3A."))

E.append(dict(fid="G2-01", status="fixed", changed_in_repair=True,
  reproduced="Re-derived in repair job 7550: base fits reproduce the released stage_extension_all_gene_results.tsv exactly (max |dlogFC| 0, all seven counts equal). F1 vs F0 without GSE130970 / GSE135251 / GSE162694: 14 / 929 / 651 (155 / 115 / 140 units); all 14 are among the 1,338 (635 and 565 for the other two). Alone: 252 / 3 / 0 (50 / 90 / 65 units). Same sign in all three separate single-cohort fits 95.4%; 95.8% from the cohort coefficients of the interaction model (job 1032, checker job 9050). F4 vs F0 without GSE130970: 4,800 calls, 4,596 of the 7,043 retained. New: F2-F4 vs F0 stay 1,934 to 9,546 under any one omission; adjacent steps are cohort-dependent (F1→F2 0 without GSE162694, 133 units; F2→F3 2 without GSE135251, 54 units; F3→F4 6,565 without GSE135251, 29 units).",
  what_changes="Fig 3E bar heights unchanged. Legend 3E, draft, RESULTS, STATUS and METHODS state that the F1 count is 14 without GSE130970, that F2–F4 counts stay 1,934 to 9,546 under any one omission, that later adjacent steps depend on single cohorts, and the F4 cohort heterogeneity. The round-1 library-complexity sentences (104 genes and others) are removed because H82 returned untestable and its reading admits no complexity-adjusted number. The candidate S3 F1 panel shows cohort omission and single-cohort fits only.",
  new_value="1,338 retained as released; 14 without GSE130970 (155 units); F2–F4 vs F0 1,934–9,546 under any one omission",
  old_value="1,338 presented as pooled three-cohort F1 count without qualification",
  panels=[FR + "figs3_f1_vs_f0_cohort_omission.pdf"],
  source_tables=[FR + "fibrosis_counts_cohort_omission.tsv", FR + "reproduction_vs_release.tsv", FR + "cohort_omission_summary.json", F1 + "fig3e_cohort_heterogeneity.tsv"],
  commands=REFIT + "; Rscript repair/F3ef/scripts/render_s3_stage_sensitivity_repair.R; pdftotext -layout and pdftoppm PNG viewed; build_text_edits.py (each old_text found once, no overlaps, no em dash or evidence tag)",
  text_edits=TE["G2-01"],
  note=H82 + " So the complexity text is dropped, not reworded; the checker's '104 remain' and 'singly and jointly' problems no longer apply. 95.8% is attributed to the cohort-by-stage model in legend and draft; RESULTS also gives 95.4% from single-cohort fits. The METHODS sentence sits after the existing GNMT/MAT1A rule so 'these comparisons' keeps its F0-reference meaning. RESULTS 1128 and STATUS 817 edits are new (conflict 19). The STATUS 540-541 edit also carries R12-04's NAS0 qualifier."))

E.append(dict(fid="G2-03", status="fixed", changed_in_repair=True,
  reproduced="Released 54 genes / 216 effects / 59 marks (30/13/11/5) re-derived (round 1 and checker). Disease-only F0 from released f0-arms r8 arm B (re-derived in repair): 18 of 30 F0→F1 marks kept (all same sign); 12 only-F0→F1 genes lose the call; reselection 45 genes / 50 marks (21/13/11/5), 42 shared, FCHO2, HPS6, SEC14L2 enter. Repair job 7550, marks kept (q<0.05, same sign) without GSE130970 / GSE135251 / GSE162694 per contrast F0→F1, F1→F2, F2→F3, F3→F4: 0/6/6/0, 15/3/1/2, 12/0/3/4; smallest F0→F1 q without GSE130970 0.103; only OMG (F2→F3) holds under all three omissions.",
  what_changes="Fig 3F membership and the installed panel (black dots) are unchanged. Legend 3F and the draft now state per-contrast FDR with the 15-gene split-budget set (merged with F3cdg R16-02, conflict 22), that 29 genes enter only through F0→F1 marks, that none of those 30 marks keeps q<0.05 without GSE130970 (18 do with disease-only F0), and that none of the 13 F1→F2 marks does without GSE162694. RESULTS and STATUS carry the counts. The filled/open candidate panel is withdrawn: its second check used joint library-complexity covariates, and H82 returned untestable.",
  new_value="54 genes / 216 effects / 59 marks unchanged; cohort-omission and disease-only-F0 qualifiers in legend and text",
  old_value="59 marks (30/13/11/5) with no cohort or F0-composition qualifier",
  panels=[],
  source_tables=[FR + "fig3f_marks_cohort_omission.tsv", FR + "fibrosis_counts_cohort_omission.tsv", B + "figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_results.tsv.gz", F1 + "fig3f_mark_sensitivity.tsv"],
  commands=REFIT + "; rnaseq python reselection from f0_arm_results.tsv.gz arm B with gene_level_coloc_tier12.csv (436 genes, login node); build_text_edits.py",
  text_edits=TE["G2-03"],
  note=H82 + " The checker's key fix ('q<0.05, not under both checks (33)') was applied and rendered before the H82 result appeared (PNG checked); the PDF is kept in repair/withdrawn_H82_untestable/ as a record only and must not be installed. The legend and fig3.md edits replace F3cdg R16-02 E019/E020; drop those two at install. Legend 3F now describes only the installed panel, so the checker's legend-panel mixing problem no longer applies to 3F. Omission removes 12% to 69% of a contrast's units, so part of the mark loss reflects lower power [Inference]; the draft says this plainly."))

r4 = dict(first["R03-04"]); r4["changed_in_repair"] = False
E.append(r4)

E.append(dict(fid="R03-05", status="fixed", changed_in_repair=True,
  reproduced=first["R03-05"]["reproduced"] + " Re-derived again in repair from the released f0_arm_results.tsv.gz: same 18/30, same 12 genes, 45 genes / 50 marks, 42 shared, 3 enter.",
  what_changes="Legend 3F and the draft state that 18 of the 30 F0→F1 marks keep q<0.05 with disease-only F0 (G2-03 edits). RESULTS names the 12 genes that lose the call and the disease-only reselection (45 genes, 50 marks; FCHO2, HPS6 and SEC14L2 enter). The filled/open candidate panel that combined this check with complexity covariates is withdrawn (H82 untestable).",
  new_value=first["R03-05"]["new_value"], old_value=first["R03-05"]["old_value"],
  panels=[],
  source_tables=[B + "figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_results.tsv.gz", F1 + "fig3f_mark_sensitivity.tsv"],
  commands="rnaseq python on released f0_arm_results.tsv.gz arm B joined to the released 3F source table (login node)",
  text_edits=[],
  note="Text is in G2-03's legend 3F, fig3.md 124-127, RESULTS 37-38 and STATUS 110 edits."))

E.append(dict(fid="R03-06", status="fixed", changed_in_repair=True,
  reproduced=first["R03-06"]["reproduced"] + " Repair: arm C n_reference = 103 and all 12 arm counts asserted in the render script; F1 omission counts from job 7550.",
  what_changes="The repaired S3 F0-reference panel prints all 12 counts under the points and gives n for all three references (n=103 added for the control-status arm). The S3 F1-versus-F0 panel shows the released fit, the three omissions and the three single-cohort fits (job 7550 values); the complexity bars are removed (H82 untestable). The legend 3E and S3 legend edits name these panels and drop complexity wording; the S3 legend edit is anchored at its actual lines 175-176.",
  new_value="S3 panels present: F0-reference arms with printed counts and n; F1 vs F0 by cohort omission and single-cohort fit",
  old_value="claimed but absent",
  panels=[FR + "figs3_f0_reference_arms.pdf", FR + "figs3_f1_vs_f0_cohort_omission.pdf"],
  source_tables=[FR + "figs3_f0_reference_arms.tsv", FR + "figs3_f1_vs_f0_cohort_omission.tsv"],
  commands="Rscript repair/F3ef/scripts/render_s3_stage_sensitivity_repair.R (adapted copy of ws/F3ef/render_s3_stage_sensitivity.R); pdftotext -layout shows 1,338/4,266/6,190/7,043, 783/3,946/6,459/6,929, 814/3,679/6,166/6,640 and all three n; pdftoppm PNGs viewed",
  text_edits=TE["R03-06"],
  note="Install both legend edits only together with repair/figs3_f0_reference_arms.pdf and repair/figs3_f1_vs_f0_cohort_omission.pdf. The round-1 F3ef/figs3_f0_reference_arms.pdf and F3ef/figs3_f1_vs_f0_cohort_complexity.pdf are superseded. H82's prespec targeted the same legend 145-146 old_text, but its result proposes no text (untestable), so there is no collision."))

E.append(dict(fid="R12-04", status="fixed", changed_in_repair=True,
  reproduced="Checker was right. Re-derived from the release: stage_extension_sample_manifest.tsv (NAS_adjacent rows) joined to five_cohort_sample_manifest.tsv group_binary, counts from stage_extension_all_gene_results.tsv. NAS0 46 = GSE130970 4 + GSE135251 10 + GSE162694 31 controls + 1 GSE162694 disease. NAS1-2 65 = 63 disease + 2 GSE130970 controls. NAS3-4 117 (32/61/24) and NAS5-8 167 (22/105/40) are all disease. NAS1-2→3-4: 3 genes (2 up, 1 down). NAS3-4→5-8: 4,305 (1,954 up, 2,351 down). Versus NAS0: 5,338 / 9,337 / 9,841.",
  what_changes="Legend S3 and RESULTS describe the versus-NAS0 counts as mainly disease-versus-control. The round-1 clause 'the only adjacent step within disease, NAS 1–2 to 3–4, has 3 genes' was false and is replaced by both within-disease steps: 3 genes (NAS 1–2→3–4, with 2 controls in NAS 1–2) and 4,305 genes (NAS 3–4→5–8).",
  new_value="NAS0 = 45 controls + 1 disease; within disease NAS 1–2→3–4 = 3 genes, NAS 3–4→5–8 = 4,305 genes",
  old_value="'continuous NAS activity is centered on NAS0 (46 units with fixed scores)' with NAS0 counts read as activity grade; round-1 draft said NAS 1–2→3–4 was the only within-disease step",
  source_tables=[F1 + "nas0_reference_composition.tsv", B + "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_sample_manifest.tsv", B + "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"],
  commands="rnaseq python merge of stage_extension_sample_manifest.tsv and five_cohort_sample_manifest.tsv plus BH counts by sign from stage_extension_all_gene_results.tsv (login node)",
  text_edits=TE["R12-04"],
  note="STATUS 540-541 NAS0 qualifier is merged into G2-01's STATUS edit. RESULTS 1104 already printed 1,954 / 2,351 for NAS 3–4→5–8."))

E.append(dict(fid="G2-04", status="fixed", changed_in_repair=True,
  reproduced="Repair job 7550 (base reproduces release) plus released f0-arms r8. GNMT F0→F1 q 0.00913 / 0.0426 / 0.0430 (all F0 / disease-only / control-status term); without GSE130970 / GSE135251 / GSE162694 0.207 / 0.0429 / 0.0338. MAT1A F0→F1 2.45e-4 / 9.5e-4 / 1.6e-3; 0.0892 / 0.00251 / 0.00631. CYP2C19 F1→F2 0.0477 (0.052 control-status term); 0.0239 / 0.219 / 0.499; CYP2C19 F1 vs F0 0.349. F2-F4 vs F0 for all three: negative in every fit; max q 2.8e-4 over F0 references, 3.8e-4 without GSE130970; without GSE135251 MAT1A F4 0.0592 (11 F4 units), CYP2C19 F2 0.0338, GNMT F4 0.0238; without GSE162694 max 0.0061. The checker was right that 'in every fit' overclaimed.",
  what_changes="Draft fig3.md, RESULTS 1129-1131 and 557-559, and PAPER 253-254 anchor the trio on F2–F4 vs F0 and report the earlier differences with their cohort dependence: GNMT and MAT1A at F1 do not hold without GSE130970; CYP2C19 F1→F2 does not hold without GSE135251 or GSE162694. Wording is cross-sectional ('lower abundance at F1 than at F0'). The round-1 sentence 'None of these early steps is a current stage call' is removed because those q values are below 0.05 in the released fit. Complexity-covariate q values are removed (H82 untestable).",
  new_value="F2–F4 vs F0 lower for all three in every fit (q<0.001 under F0 references and without GSE130970; ≤0.034 otherwise except MAT1A F4 without GSE135251, q 0.059); F1 differences of GNMT/MAT1A depend on GSE130970",
  old_value="GNMT decreases from F0 to F1 (q=0.00913 / 0.0426 / 0.0430); MAT1A F0 to F1; CYP2C19 F1 to F2 (q=0.0477); RESULTS 557 'stays supported'",
  source_tables=[FR + "hero_genes_cohort_omission.tsv", B + "figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_hero_genes.tsv"],
  commands=REFIT + "; rnaseq python pivot of f0_arm_results.tsv.gz for GNMT/MAT1A/CYP2C19; build_text_edits.py",
  text_edits=TE["G2-04"],
  note="The METHODS rule ('Early-stage GNMT and MAT1A claims that do not survive these comparisons are not retained') refers to the F0-reference comparisons, which GNMT and MAT1A pass, so the text reports cohort dependence instead of invoking it. PAPER now says '(the GNMT and MAT1A differences at F1 depend on GSE130970)' per the checker; CYP2C19 has no F1 difference (q 0.349). Intro 93-94 and fig5_discussion 173 ('cross-sectional fibrosis-stage abundance') stay true for F2–F4 vs F0. RESULTS 557-559 are the conflict-19 edits. The round-1 complexity values (GNMT 0.360, CYP2C19 0.0729 joint) are review-only after H82."))

E.append(dict(fid="G2-05", status="covered_by_plan", changed_in_repair=True,
  reproduced="Values reproduced in round 1 and by the checker: F3→F4 281 released, 643 log2 library size, 115 Gini, 117 norm factor, 166 joint; H82 D1 gives the same counts. H82 is the plan for gate G2-02 (complexity covariates).",
  what_changes="No change. H82 returned UNTESTABLE (S5: PC-3 failed), and its reading puts no complexity-adjusted count in the paper, so RESULTS 1123 keeps 281 and the round-1 range edit is withdrawn. The cohort-omission F3→F4 counts (102 / 6,565 / 508 without GSE130970 / GSE135251 / GSE162694; 82 / 29 / 77 units) are recorded in G2-01's RESULTS 1112 edit.",
  new_value="281 (unchanged)", old_value="281",
  source_tables=[B + "Analysis/Hypotheses/H82-20261015T125052Z/RESULT_H82.md", FR + "fibrosis_counts_cohort_omission.tsv"],
  commands="read RESULT_H82.md (verdict computed by the H82 run, job 6220)",
  text_edits=[],
  note=H82))

r7 = dict(first["R03-07"]); r7["changed_in_repair"] = True
r7["blocked_on"] = "H51 result (RESULT_H51.md; stage 2 job 22252470 still running at 17:33Z). The job-1032 disclosure is now recorded in Analysis/Hypotheses/H51-20261015T112617Z/DISCLOSURE_ADDENDUM_H51.md."
r7["what_changes"] = "No change in this round. Plan H51 tests replication of the 3E/3F stage effects in the GSE213621 disease units; draft wording waits for RESULT_H51.md. The disclosure that job 1032 had already fitted and reported the GSE213621 F3-F4 vs F0-F1 contrast is now in DISCLOSURE_ADDENDUM_H51.md, which says H51 cannot be presented as a prediction confirmed on unseen data."
r7["note"] = first["R03-07"]["note"] + " Repair: H51 status checked in squeue and the run directory; no RESULT_H51.md yet."
E.append(r7)

json.dump(E, open(RD + "repair/F3ef/entries_repair.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(E, ensure_ascii=False))
