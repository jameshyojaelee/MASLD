"""Assemble the F4abcd StructuredOutput payload from edits.json plus per-finding records."""
import json
from pathlib import Path

HERE = Path(__file__).parent
WS = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F4abcd"
FIG = "figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd"
H59 = "Analysis/Hypotheses/H59-20261015T112432Z"
edits = json.loads((HERE / "edits.json").read_text())

P4C = [f"{FIG}/fig4c_hotspot_stage_heatmap_H59_candidate.pdf",
       f"{FIG}/figs4c_hotspot_stage_heatmap_full_H59_candidate.pdf"]
ST4C = [f"{FIG}/source_tables/fig4c_h59_main_24rows.tsv",
        f"{FIG}/source_tables/fig4c_h59_registry_117.tsv",
        f"{FIG}/source_tables/fig4c_h59_registry_117_with_familywise_interval.tsv",
        f"{H59}/out/outcome/results_primary_family.tsv"]
CMD4C = (f"login node, saved tables only: {WS}/rederive/rederive.py (log rederive.log); "
         f"{FIG}/producers/01_build_fig4c_h59_table.py; Rscript {FIG}/producers/02_render_fig4c_h59.R; "
         "pdftotext -layout on installed fig4c_hotspot_stage_heatmap.pdf and the candidate; pdftoppm PNG check in png_check/")

meta = {
 "R04-01": dict(
  status="fixed",
  reproduced=("62-donor OLS refit reproduces the installed call exactly: Ductular injury beta 0.402086, p 0.000322, "
              "q 0.037633, HC3 q 0.155826, 1/117 calls; installed Fig 4C pdftotext shows the only '*' on '20 · Ductular injury'. "
              "Dropping the two GSE189600 donors gives q 0.098523 and 0/117 calls. H59 (GLS varIdent(dataset), 60 donors, 19,999 "
              "within-dataset permutations, seed 20260923; gates G0-G6 pass) gives 0/117 selected: Ductular injury beta 0.1365 "
              "(0.0121 to 0.2609), perm p 0.0367, q 0.323375; Stromal ECM 0.0837 (-0.0672 to 0.2346), q 0.577711."),
  what_changes=("Fig 4C loses its only BH star and the Ductular injury beta drops from 0.402 to 0.137. Abstract, 4C legend, "
                "PAPER claim 2, RESULTS, STATUS, README, ROADMAP, intro and fig5_discussion change from '1/117 primary-selected "
                "(Ductular injury q=0.03763)' to '0/117 selected (smallest q 0.323375)'. The 62-donor q 0.037633 becomes a "
                "supplementary sensitivity that depends on equal variances and on two GSE189600 donors without ambient correction. "
                "Stromal ECM text adds 'cannot exclude the printed effect' (PC4 detection 0.794)."),
  old_value="1/117 primary-selected; Ductular injury beta 0.4021, q 0.03763, HC3 q 0.1558 (62 donors, OLS)",
  new_value=("0/117 selected; Ductular injury GLS beta 0.1365 (95% CI 0.0121 to 0.2609), perm q 0.323375; Stromal ECM 0.0837 "
             "(-0.0672 to 0.2346), q 0.577711 (60-donor universe, fitted donors in GSE202379 and GSE244832)"),
  panels=P4C, source_tables=ST4C + [f"{FIG}/source_tables/R04-01_02_lodo_family.tsv"], commands=CMD4C,
  note=("Text edits here also carry the R04-02 slopes, R04-03 donor counts and G1-01 gate sensitivity. Under the H59 E9 row "
        "rule Stromal ECM (hepatocytes 8) drops out of the 24 main 4C rows; it stays in S4C, and fig4.md says so. "
        "PAPER:336-337 and fig5.md:41-46 ('robust', 'equally robust') belong to R06-02 and are not drafted here. Gene Catalog r8 "
        "numbers do not change (H59 destination section); the passport_program_context basis becomes 'prespecified programs'. "
        "The text treats H59 as the primary fit, which needs user adoption. Nothing is installed.")),
 "R04-02": dict(
  status="fixed",
  reproduced=("Within-dataset slopes on the 62-donor design: Ductular injury GSE202379 0.1074 (SE 0.0633, p 0.0986, n 37), "
              "GSE244832 0.7423 (SE 0.2568, p 0.0106, n 18); Stromal ECM 0.0373 (p 0.633) and 0.7556 (p 0.0089). Omit GSE244832: "
              "0 calls, Ductular beta 0.166314, q 0.425144. Omit GSE202379: 0 calls, beta 0.731994, q 0.172025. H59 Cochran Q p "
              "0.0164 (Ductular) and 0.00676 (ECM); GLS CI coverage 0.942 and 0.955 at true beta 0 and 0.30 (G3 pass)."),
  what_changes=("The printed 4C beta for Ductular injury goes from 0.402 to 0.137 and for Stromal ECM from 0.352 to 0.084 once "
                "each dataset gets its own residual variance. METHODS now states the GLS variance model and permutation test, "
                "and the text reports both dataset slopes and says the pooled beta is a precision-weighted average."),
  old_value="Ductular injury beta 0.402 (OLS, equal variance); Stromal ECM 0.352",
  new_value="Ductular injury 0.1365; Stromal ECM 0.0837 (GLS, per-dataset variance)",
  panels=P4C, source_tables=[f"{FIG}/source_tables/R04-01_02_lodo_family.tsv", f"{FIG}/source_tables/S4_native_replication_and_omission_corrected.tsv", f"{H59}/out/outcome/named_programs.tsv"],
  commands=f"{WS}/rederive/rederive.py (login node, saved tables)",
  note=("The fig4.md and RESULTS text carrying the slopes is in the R04-01 entry. My omit-GSE202379 full-family refit gives "
        "Ductular q 0.172025, where the G1-01 verifier quoted 0.118 [Inference: theirs was a single-program refit]; the "
        "conclusion is the same.")),
 "R04-03": dict(
  status="fixed",
  reproduced=("Exact match. Donors in datasets with more than one stage, 62-donor design: hepatocytes 57/62 in 3 datasets "
              "(GSE189600, GSE202379, GSE244832); cholangiocytes 38/43, fibroblasts 43/47, macrophages 25/30, T cells 53/58, each "
              "in 2 datasets. GSE174748 (2 Healthy) and GSE185477 (3 Healthy) are single-stage. Macrophage Steatosis donors 4 (1+3). "
              "H59 fitted n equals these counts, with hepatocytes 55 after GSE189600 is dropped."),
  what_changes=("'62 biological donors across five datasets' is replaced, for the stage slope, by 57 hepatocyte donors in three "
                "datasets and 38/43/25/53 donors in two datasets (62-donor fit), and by 55/38/43/25/53 donors in GSE202379 and "
                "GSE244832 under H59 (fig4.md, RESULTS, 4C legend)."),
  old_value="62 donors across five datasets", new_value="stage-informative donors: 55 hepatocyte, 38 cholangiocyte, 43 fibroblast, 25 macrophage, 53 T-cell, all in 2 datasets (H59)",
  source_tables=[f"{FIG}/source_tables/R04-03_informative_donors.tsv", f"{WS}/rederive/R04-03_design_by_lineage_dataset_stage.tsv"],
  commands=f"{WS}/rederive/rederive.py",
  note="The fig4.md and 4C legend edits that print the H59 counts sit in the R04-01 and R04-05 entries."),
 "R04-04": dict(
  status="partly_reproduced",
  reproduced=("All numbers reproduce: ECM raw 0.350898 to corrected 0.127891 (-63.6%), delta -0.223006, nominal p 0.006417 "
              "(HC3 0.024386), 117-program BH q 0.478071; two-program BH q 0.0128 (HC3 0.0488); ECM delta p is the smallest of 106 "
              "computable programs. Ductular delta +0.016422 (p 0.671). Not reproduced as a defect: METHODS:666-668 prespecifies "
              "the 117-program family for paired changes, so calling it in a two-program family would be post hoc."),
  what_changes=("No call changes; the ECM paired change stays unsupported at the prespecified 117-program BH. The text now prints "
                "the 64% reduction and nominal p 0.0064 (rank 1 of 106) beside q 0.4781, plus the H59 E4 ambient-corrected GLS betas "
                "(ECM 0.0518, q 0.559; Ductular 0.1428, q 0.270). With no program selected under H59, the existing rule "
                "(common.py:235-236: not primary-selected gives not_applicable_not_frozen_disease_association) makes the "
                "attenuation classification inapplicable to Ductular injury as well."),
  old_value="ECM delta BH q 0.4781 printed alone", new_value="ECM -63.6%, nominal p 0.0064 (HC3 0.024), smallest of 106; BH q 0.4781",
  source_tables=[f"{FIG}/source_tables/R04-04_focal_ambient_rows.tsv", "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/ambient_program_effects.tsv", f"{H59}/out/outcome/named_programs.tsv"],
  commands=f"{WS}/rederive/rederive.py",
  note="The Ductular injury classification change is a mechanical application of the existing rule to the H59 result [Inference; the ambient producer was not rerun]."),
 "R04-05": dict(
  status="fixed",
  reproduced=("Installed 4C draws beta ± 1.96 × OLS SE: 11/24 intervals exclude 0, with 1 BH call and 0 HC3 calls. Under H59, "
              "17/24 nominal GLS intervals exclude 0 (19/117 overall) while 0/117 of the 117-program intervals (beta ± per-program "
              "permutation critical |t| at 0.05/117 × SE) exclude 0. Healthy-column shift when single-stage datasets are dropped: "
              "median |shift| 0.070, max 0.714; median |Steatohepatitis-Healthy| contrast 0.272 installed versus 0.410 on the two "
              "multi-stage datasets."),
  what_changes=("Fig 4C shows GLS betas with a black nominal 95% CI and a gray 117-program interval, no star, and heatmap columns "
                "built only from GSE202379 and GSE244832. A reader sees 0/24 gray intervals excluding zero, matching 0/117 family "
                "calls. The legend states what each bar means and the fitted donor counts."),
  old_value="11/24 nominal OLS CIs exclude 0; 1 BH star; Healthy column includes single-stage datasets centered within themselves",
  new_value="17/24 nominal GLS CIs exclude 0; 0/24 117-program intervals exclude 0; no star; heatmap from the two multi-stage datasets",
  panels=P4C, source_tables=ST4C + [f"{FIG}/source_tables/R04-05_healthy_column_single_stage_datasets.tsv", f"{H59}/out/null_critical_abs_t.tsv"],
  commands=CMD4C,
  note=("The heatmap departs from H59 E9's stage means on 60 donors: it uses fit donors only, which is the R04-05 fix. The 24 "
        "rows equal H59 E9's row choice (asserted in 01_build). The legend edit also carries the G1-03 T-cell disclosure. "
        "Font 6, no bold, useDingbats=FALSE, zero line #9E9E9E; PNG checked.")),
 "R04-06": dict(
  status="fixed",
  reproduced=("Installed S4_native_replication_and_omission.pdf prints the 64-donor design: a 'GSE174748 (n=4)' row and omission n "
              "46, 27, 62, 61, 60. On the corrected design GSE174748 and GSE185477 have fewer than two stage levels and GSE189600 has "
              "nonpositive residual df, so they have no dataset row; omission n are 60, 59, 60, 25, 44."),
  what_changes=("S4 native forest now has 8 rows: pooled n=62 (ECM 0.3515, Ductular 0.4021), GSE202379 n=37 (0.0373, 0.1074), "
                "GSE244832 n=18 (0.7556, 0.7423), and omissions at n 60/59/60/25/44 (omit GSE202379: 0.7034, 0.7320; omit GSE244832: "
                "0.0715, 0.1663), with HC3 95% CIs. The GSE174748 row is removed."),
  old_value="GSE174748 (n=4) row; omission n 46/27/62/61/60", new_value="no GSE174748 row; omission n 60/59/60/25/44",
  panels=[f"{FIG}/S4_native_replication_and_omission_corrected_candidate.pdf"],
  source_tables=[f"{FIG}/source_tables/S4_native_replication_and_omission_corrected.tsv",
                 "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/context_effects.tsv",
                 "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv"],
  commands=f"python {FIG}/producers/03_render_s4_native_replication.py (log 03_render.log); pdftotext -layout; PNG check",
  note="The panel is labelled as the 62-donor equal-variance sensitivity, consistent with H59 moving that fit to the supplement."),
 "R04-07": dict(
  status="fixed",
  reproduced=("V3 v3_primary_family.tsv: Ductular injury C1 0.694 (0.062 to 1.326), p 0.0347, Holm 0.139 (n 13); C2 1.080 "
              "(-0.091 to 2.250), p 0.0677, Holm 0.203 (n 17); Stromal ECM C1 0.300 (-0.709 to 1.309), C2 0.176 (-0.840 to 1.192), "
              "Holm 1. Model score ~ stage + sex + age (C1) or disease + sex + age (C2), HC3, Holm over 4. grep finds GSE256398 only "
              "at STATUS:323, not in the manuscript, RESULTS or PAPER."),
  what_changes=("Fig 4 text, RESULTS and METHODS gain the only new-donor test of the two programs: no test passes Holm, power is "
                "low at 13 and 17 donors, and no positive control was run. No panel number changes."),
  old_value="not reported", new_value="V3: 0/4 tests pass Holm; Ductular injury stage Holm p 0.139",
  source_tables=["Analysis/Hypotheses/V3-20261008T015239Z/results/v3_primary_family.tsv"],
  commands="column -t v3_primary_family.tsv; grep -rn GSE256398 docs/",
  note=("V3 has no RESULT file. Adding GSE256398 also needs a data-availability and dataset-roster line, not drafted here. "
        "GSE289173 was not opened.")),
 "R04-08": dict(
  status="fixed",
  reproduced=("Installed 4D (pdftotext) is the discovery-by-content matrix (hepatocyte row 20, 2, 4, 1, 2, 1; 80 on the diagonal, "
              "37 off). The ten-contrast localization result (RESULTS:679-686: ECM in fibroblasts beta 2.271, BH q 3.38e-26, 53 "
              "donors; Ductular in cholangiocytes 1.574, q 2.04e-23, 51 donors) is not on 4D; the legend (243-247) already places it "
              "in the supplement."),
  what_changes="fig4.md:107-111 now describes the 80/37 matrix shown in 4D and cites the localization contrasts to Figure S4 with their betas and donor counts.",
  old_value="text says 4D shows ten donor-paired contrasts", new_value="text says 4D shows 80/37 discovery-versus-content placement",
  source_tables=["figures/main/fig4_singlecell_programs/CANDIDATE_PANEL_INDEX.tsv"],
  commands="pdftotext -layout on the installed 4D panel; sed -n 679,700p docs/RESULTS.md",
  note="RESULTS:695-697 still records the 2026-08-17 approval of localization as 4D; that is history and is not edited."),
 "R04-09": dict(
  status="fixed",
  reproduced=("Full content call agrees with the published call for 117/117, 80 matches. Within-lineage bootstrap, seed 20261015, "
              "200 draws reproduces the finding: 75 stable matches, 34 stable mismatches, 8 unstable, per-draw 75 to 80. Joint donor "
              "bootstrap (donor resampled across all lineages), seed 20260923, 1,000 draws: 76/34/7, per-draw median 79 (central 95% "
              "77 to 80, range 74 to 82). Within-lineage, seed 20260923, 1,000 draws: 75/32/10. Both focal programs keep their call "
              "in 100% of draws in all runs. Stable means the same call in at least 95% of draws."),
  what_changes="The 80/37 counts stay as point estimates; fig4.md and the 4D legend add 76 stable matches, 34 stable mismatches and 7 unstable (10 under within-lineage resampling). Focal-program conclusions do not change.",
  old_value="80/37 with no uncertainty", new_value="80/37; 76 stable, 34 stable mismatch, 7 unstable (joint donor bootstrap, B=1,000)",
  source_tables=[f"{FIG}/source_tables/R04-09_content_lineage_bootstrap_summary.tsv", f"{FIG}/source_tables/R04-09_content_lineage_bootstrap_joint_seed20260923_B1000.tsv"],
  commands=f"sbatch job 7061 (slurm 22253009, cpu, nslab, 2 CPU, 16G): {WS}/jobs/7061/content_bootstrap.py; COMPLETED 00:15:40; pip_freeze.txt written",
  note="The joint scheme is the primary one because the donor, not the donor-lineage pair, is the unit."),
 "R04-10": dict(
  status="fixed",
  reproduced=("stromal_ecm_ductular_by_arm.tsv, arm new774_plus_GSE174478, Stromal ECM: fibrosis unadjusted beta 0.3255, q "
              "0.0452; composition-adjusted 0.1461, q 0.2303; NAS q 0.6077 unadjusted and 0.6980 adjusted (BH over 224). "
              "Single-cell: 62-donor q 0.085458; H59 q 0.577711."),
  what_changes="Intro:97-99 no longer calls ECM a 'disease-associated program'; it states that the ECM disease association rests on unadjusted bulk fibrosis (q 0.045) and is lost after composition adjustment (q 0.230).",
  old_value="'a disease-associated program'", new_value="bulk fibrosis q 0.045 unadjusted, 0.230 composition-adjusted",
  source_tables=["figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/stromal_ecm_ductular_by_arm.tsv"],
  commands="grep new774_plus_GSE174478 stromal_ecm_ductular_by_arm.tsv"),
 "R10-01": dict(
  status="fixed",
  reproduced=("Independent rerun (job 7060, the review script copied unchanged, script seeds 20260812 and 20260812+100000, 200 "
              "draws, floor 1/201 = 0.004975) matches the review exactly. Arm A, 469 records: fibrosis 26 vs null mean 18.445 "
              "(5.975 to 31.075), p 0.1343; NAS 4 vs 4.225, p 0.5124 (equals the printed values). Arm C, 457 records (the Fig 4E "
              "arm): fibrosis 20 vs 26.91 (13.975 to 38.025), p 0.8706; NAS 2 vs 2.565, p 0.6866. Arm B, 364 units: 1 vs 1.24, p "
              "0.4726; 0 vs 0.19, p 1.0. BH family 228 in every arm."),
  what_changes=("fig4.md:193-195 changes from 'fibrosis 26/114 versus null mean 18.445, P=0.134' to '20/114 versus 26.9, P=0.871' "
                "on the 457-record Fig 4E arm; RESULTS labels 26/114 as the 469-record roster. The conclusion (no exceptional "
                "specificity) holds and is stronger: Hotspot programs get fewer fibrosis calls than matched random sets."),
  old_value="26/114 vs 18.445, P=0.134; NAS 4 vs 4.225, P=0.512", new_value="20/114 vs 26.91, P=0.871; NAS 2 vs 2.565, P=0.687 (457 records)",
  source_tables=[f"{FIG}/source_tables/R10-01_vocabulary_null_summary_all_arms.tsv", f"{WS}/jobs/7060/out/summary_all_arms.tsv"],
  commands=f"sbatch job 7060 (slurm 22253000, cpu, nslab, 4 CPU, 16G): Rscript {WS}/jobs/7060/t3_hotspot_vocab_null_774.R; COMPLETED 00:06:05; sessionInfo.txt written",
  note="A weighted-score null was not run. HIGH_IMPACT_PATTERNS.md:76 also prints 26/114 vs 18.4 but labels it 469 donors, so it is left alone."),
 "R13-04": dict(
  status="partly_reproduced",
  reproduced=("Reproduced: 1/117 primary and 0/117 HC3 on the 62-donor fit, 0/117 under H59, 37/117 lineage mismatch; abstract:48-49 "
              "and intro:55-56 say the registry 'resolved' multicellular remodeling. Not reproduced: the finding's '0/224 in bulk Fig "
              "4E'. The installed legend and fig4.md:116-118 print 16 fibrosis and 1 NAS association at BH q<0.05 (17/224), and 3 "
              "fibrosis after composition adjustment. The 26/114 vocabulary number is replaced by R10-01."),
  what_changes=("Abstract (edit in R04-01), intro:54-56, fig1.md:44-45 and PAPER idea 2 (122-123) stop saying the single-cell "
                "programs 'resolve' multicellular remodeling and say instead that they score fixed programs in five lineages, with no "
                "family-level stage association. Title clause 2 is not decided: six conditional drafts follow the H52/H54 joint rule, "
                "two 'kept' rows (covariance of single-cell-derived programs in bulk liver, per fixed condition 1) and four 'removed' "
                "rows (the review's 'maps cross-sectional multicellular remodeling')."),
  old_value="'resolved remodeling into multicellular programs'; title 'the coordination of multicellular remodeling'",
  new_value="'scored a prespecified 117-program registry in five cell lineages'; title clause pending H52/H54",
  blocked_on="Title clause only: H52 and H54 have no RESULT yet; the joint rule decides which conditional title draft applies.",
  source_tables=["figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv", f"{H59}/out/outcome/results_primary_family.tsv"],
  commands="grep -n 'multicellular' docs/manuscript/draft/*.md docs/PAPER.md",
  note=("The same title replacement applies to PAPER.md:64-65, README.md:102 and STATUS.md:237, which repeat the clause. The "
        "supports/untestable row keeps the same title wording as supports/supports; its one-cohort qualifier goes in the abstract "
        "and Results. The supports/fails row also needs Results text reporting a discovery-cohort result with both replications, "
        "not drafted until H52/H54 exist.")),
 "G1-01": dict(
  status="partly_reproduced",
  reproduced=("Job 1711 arms reproduce from its log: no gate q 0.037633 (1/117); 50-cell lineage gate q 0.076511 (0/117); 20-cell "
              "gate 0.037633; v2 QC plus 50-cell gate 0.061609; drop GSE189600 0.098523. Job 7062 confirms the GSE189600 donors carry "
              "43 and 99 hepatocyte nuclei (62-roster hepatocyte minimum 43). Not reproduced as stated: METHODS:27-28 makes the "
              "50-cell gate primary only for the two targeted content-lineage tests, so the 117-program family did not break its own "
              "rule; the gate result is a sensitivity."),
  what_changes=("Under H59 the call is gone at any gate. METHODS states that the 50-cell gate covers only the targeted tests and "
                "that a 50-cell gate removes the 62-donor call (q 0.0376 to 0.0765); RESULTS (R04-01 edit) reports this sensitivity."),
  old_value="gate scope unstated", new_value="117-program models use no per-donor cell minimum; 50-cell gate sensitivity q 0.076511, 0/117",
  source_tables=["Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1711/arm_focal.tsv", f"{FIG}/source_tables/G1-03_cells_per_donor_lineage_summary.tsv"],
  commands="cat jobs/1711/qc_gate_refit.log; job 7062"),
 "G1-02": dict(
  status="fixed",
  reproduced=("ena_runinfo.tsv lists 12 snRNA-seq runs, two each for healthy_rep1-3 and NASH_rep1-3; lib_donor_collapse.R:24-31 "
              "ingests one run each from 587/589/593/595 and leaves open whether replicates are distinct people; QC_HARMONIZATION "
              "provenance 'aggressive_post_QC_loss_99.5pct'; 142 nuclei in scored lineages, 131 absent from v2 QC (pass 0.077). "
              "SINGLE_CELL_HISTORY.md:40 says Percoll/EpCAM hepatocyte isolation while ENA titles say snRNA-seq [Unverified which "
              "is correct; no GEO soft file on disk]."),
  what_changes=("Under H59 GSE189600 is outside the primary universe, so no primary number depends on it. METHODS discloses the "
                "partial ingest, the 142 nuclei and the unconfirmed donor identity, and says the dataset enters only the 62-donor "
                "sensitivity, whose Ductular call (q 0.0376) disappears without it (q 0.0985)."),
  old_value="undisclosed: 2 of 6 replicate samples, one run each", new_value="disclosed in METHODS; dataset only in the 62-donor sensitivity",
  source_tables=["data/GSE189600/metadata/ena_runinfo.tsv", "Analysis/SingleCell/results_gpu_v2_phase05/QC_HARMONIZATION.tsv", "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1711/v1_cells_passing_v2_by_dataset.tsv"],
  commands="grep SRR226775 ena_runinfo.tsv; sed -n 15,40p lib_donor_collapse.R",
  note="A full re-ingest of all 12 runs was not run; under H59 it cannot change a primary number."),
 "G1-03": dict(
  status="fixed",
  reproduced=("Job 7062 reproduces exactly. 62 roster: T cells 17/58 below 20 and 28/58 below 50 (min 1); cholangiocytes 6/43 below "
              "50 (min 30); fibroblasts 10/47 (min 32); macrophages 5/30 (min 39); hepatocytes 1/62 (min 43). H59 fit: hepatocytes "
              "0/55 below 50 (min 101); T cells 17/53 below 20 and 28/53 below 50 (min 1); other lineages unchanged. Lowest T-cell "
              "H59 perm p 0.1947 (q 0.520)."),
  what_changes="No call changes in any arm. RESULTS gains a row with cells-per-donor counts under the H59 fit and the 4C legend states 17 of 53 T-cell donors have fewer than 20 T cells.",
  old_value="not reported", new_value="T cells 28/53 fitted donors below 50 cells (17 below 20, min 1)",
  source_tables=[f"{FIG}/source_tables/G1-03_cells_per_donor_lineage_summary.tsv", f"{WS}/jobs/7062/out/cells_per_donor_lineage_62roster.tsv"],
  commands=f"sbatch job 7062 (slurm 22253013, cpu, nslab, 2 CPU, 24G): {WS}/jobs/7062/cell_counts.py; COMPLETED 00:00:23"),
 "G1-04": dict(
  status="fixed",
  reproduced=("QC_HARMONIZATION.tsv v1 columns: no MT filter, min_counts 200, min_genes unknown, no doublet step, no ambient "
              "correction, for every dataset; v2: MT<20%, min_genes 200, min_counts 500, Scrublet 0.25. Job 1711 log: v2 total "
              "1,030,882 cells; pass fractions GSE185477 0.494, GSE244832 0.793, GSE174748 0.891, GSE202379 0.970, GSE189600 0.077; "
              "Spearman 0.995 between all-cell and QC-cell Ductular donor scores. METHODS single-cell section has no per-cell QC "
              "statement (grep hit only at line 94, bulk)."),
  what_changes="METHODS gains the per-cell filters of the scored atlas and the QC-sensitivity numbers. No call changes.",
  old_value="per-cell QC not stated", new_value="200-count floor, no MT filter, no doublet removal; harmonized-QC pass 7.7% to 97.0% by dataset",
  source_tables=["Analysis/SingleCell/results_gpu_v2_phase05/QC_HARMONIZATION.tsv", "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1711/v1_cells_passing_v2_by_dataset.tsv"],
  commands="cut -f1-18 QC_HARMONIZATION.tsv; grep total jobs/1711/qc_gate_refit.log"),
}

order = ["R04-01", "R04-02", "R04-03", "R04-04", "R04-05", "R04-06", "R04-07", "R04-08", "R04-09", "R04-10",
         "R10-01", "R13-04", "G1-01", "G1-02", "G1-03", "G1-04"]
entries = []
for fid in order:
    m = dict(meta[fid])
    ent = dict(fid=fid)
    ent.update(m)
    ent["text_edits"] = [dict(file=e["file"], line=e["line"], old_text=e["old_text"], new_text=e["new_text"],
                              claim_level=e["claim_level"]) for e in edits if e["fid"] == fid]
    entries.append(ent)

cross = [
 "Ductular injury q 0.037633/0.03763 (1/117 primary-selected): abstract.md:50-51, fig4.md:75-77, 05_figure_legends.md:235-237, PAPER.md:266-269, RESULTS.md:633-634, 1600, 1641, 2785-2786, STATUS.md:23 and 822, README.md:13, ROADMAP.md:19-20, manuscript/README.md:28, 01_intro.md:99-100, fig5_discussion.md:44-46; installed Fig 4C star. All replaced by H59 0/117 (smallest q 0.323375).",
 "Ductular injury beta 0.402086/0.4021 and Stromal ECM 0.351525/0.3515 as primary effects: fig4.md:75-76, legend:236, PAPER:268, RESULTS:632-633; installed Fig 4C forest. H59: 0.1365 and 0.0837.",
 "'62 donors across five datasets': fig4.md:74, METHODS.md:651-652, PAPER.md:266, RESULTS.md:627-628 and 1606-1608; stage-informative counts 55/38/43/25/53 in 2 datasets under H59.",
 "Vocabulary null 26/114 vs 18.445, P=0.134: fig4.md:194-195, RESULTS.md:2146-2147 (edited), RESULTS.md:1634-1636 (labelled), HIGH_IMPACT_PATTERNS.md:76 (already labelled 469 donors, left).",
 "80/37 lineage match: fig4.md:54-55, legend:245-246 (bootstrap added; counts unchanged).",
 "Title clause 'the coordination of multicellular remodeling': title.md:10-11, PAPER.md:64-65, README.md:102, STATUS.md:237 (conditional on H52/H54).",
 "S4 native forest omission n 46/27/62/61/60 and GSE174748 (n=4) row: installed S4_native_replication_and_omission.pdf only; candidate prints 60/59/60/25/44.",
]
payload = dict(ws="F4abcd", entries=entries, cross_figure_numbers=cross)
(HERE / "structured_output.json").write_text(json.dumps(payload, indent=1, ensure_ascii=False))
print(len(entries), sum(len(e["text_edits"]) for e in entries))
