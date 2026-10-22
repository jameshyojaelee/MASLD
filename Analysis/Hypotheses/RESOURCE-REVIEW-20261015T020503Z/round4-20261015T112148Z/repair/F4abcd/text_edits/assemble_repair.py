"""Assemble the F4abcd repair payload: first-round entries + repaired edits.json + repair notes."""
import json
from pathlib import Path

HERE = Path(__file__).parent
R4 = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z")
FIG = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd"
REP = R4 / "repair/F4abcd"
first = {e["fid"]: e for e in json.load(open(R4 / "repair_inputs/F4abcd.json"))["entries"]}
edits = json.load(open(HERE / "edits.json"))
by_fid = {}
for e in edits:
    by_fid.setdefault(e["fid"], []).append({k: e[k] for k in ("file", "line", "old_text", "new_text", "claim_level")})

P4C_REP = [f"{FIG}/repair/fig4c_hotspot_stage_heatmap_H59_candidate_repair.pdf",
           f"{FIG}/repair/figs4c_hotspot_stage_heatmap_full_H59_candidate_repair.pdf"]
REND = (f"sbatch jobs 7560 (slurm 22256232) and 7561 (slurm 22256270), cpu, nslab, 1 CPU, 4G: Rscript "
        f"{FIG}/repair/producers/02r_render_fig4c_h59.R (log 02r_render.log, sessionInfo_02r_render.txt); "
        "pdftotext -layout and pdftoppm PNG in repair/png_check/. Row labels identical to the first-round candidates "
        "(24 and 117 rows, no star); familywise table byte-identical after sort.")
INSTALL = ("The text edits and the repaired Fig 4C/S4C PDFs must be installed in the same step, so no reader sees H59 "
           "text beside the installed 62-donor panels. PAPER and RESULTS now both cite H59 without a 'candidate' label, "
           "because the H59 verdict is final under plan_results.json; final adoption still follows the ROADMAP checklist "
           "and user approval.")

upd = {}
upd["R04-01"] = dict(changed_in_repair=True, panels=P4C_REP,
    commands=first["R04-01"]["commands"] + "; repair: " + REND + " Edit check: repair/F4abcd/text_edits/edits.py (41 edits, each old_text once in the current file and still once after the earlier edits to the same file are applied).",
    reproduced=first["R04-01"]["reproduced"] + " Repair: H59 results_primary_family.tsv also gives the 62-donor OLS fit a within-dataset permutation q of 0.12285 for Ductular injury (0/117 calls).",
    what_changes=first["R04-01"]["what_changes"] + " Repair: RESULTS and fig4.md add that the 62-donor call also fails the within-dataset permutation null (q 0.123); the RESULTS H59 link loses its 'candidate' label to match PAPER; STATUS says the repaired panels replace the 62-donor panels with this text; the Fig 4 README 4C row (conflict 12) changes from 'Ductular injury is the only selected primary call' to no program selected under H59.",
    note=("Text edits here also carry the R04-02 slopes, R04-03 donor counts and G1-01 gate sensitivity. Under the H59 E9 row "
          "rule Stromal ECM (hepatocytes 8) drops out of the 24 main 4C rows; it stays in S4C. " + INSTALL +
          " The RESULTS:1767 edit (E293) is left as drafted; conflict 23 (owned elsewhere) replaces it with F5ef E068. "
          "Shared line: 01_intro:99 holds the end of the R04-10 old_text and the start of this entry's old_text; the two "
          "old_texts do not overlap, so apply both as string replacements (verified sequentially). PAPER:336-337 and "
          "fig5.md:41-46 belong to R06-02. Gene Catalog r8 numbers do not change."),
    source_tables=first["R04-01"]["source_tables"] + [f"{FIG}/repair/source_tables/fig4c_h59_registry_117_with_familywise_interval.tsv"])
upd["R04-02"] = dict(changed_in_repair=True, panels=P4C_REP,
    reproduced=first["R04-02"]["reproduced"] + " Repair: H59 D2 diagnostics give GSE244832/GSE202379 residual SD ratios of 3.599 (Ductular injury) and 2.763 (Stromal ECM), variance ratios 12.95 and 7.64; across 117 programs the SD ratio has median 1.600 and range 0.067 to 15.33. PRESPEC_H59.md:163-165 defines p = (1 + count)/(1 + successful draws) = /20,000; every perm_p x 20,000 is an integer and no draw failed (G5).",
    what_changes=first["R04-02"]["what_changes"] + " Repair: METHODS now says the residual standard deviation ratio is 3.6 and 2.8 (variance ratios about 13 and 8; median SD ratio 1.6, range 0.07 to 15.3 over 117), not 'variance differed about threefold', and gives the p-value as (1 + count)/20,000.",
    commands=first["R04-02"]["commands"] + "; repair: pandas on H59 out/D2_D3_D4_D5_D6_diagnostics.tsv, D2_sd_ratio_summary.tsv and out/outcome/results_primary_family.tsv; grep on PRESPEC_H59.md",
    note=first["R04-02"]["note"])
upd["R04-03"] = dict(changed_in_repair=True,
    what_changes=first["R04-03"]["what_changes"] + " Repair: 'the slope rests on' becomes 'the slope is estimated from'.")
upd["R04-04"] = dict(changed_in_repair=True,
    new_value="ECM -63.6%, nominal p 0.0064 (smallest ordinary p of 106; HC3 p 0.024, eighth by HC3 p); BH q 0.4781",
    reproduced=first["R04-04"]["reproduced"] + " Repair: in the complete-case common universe (106 computable), ECM is first by ordinary delta p and eighth by HC3 delta p; cholangiocyte module 21, Inflammatory / oxidative (SOD2), is first by HC3 p (0.00307).",
    what_changes=first["R04-04"]["what_changes"] + " Repair: 'the smallest of 106' becomes 'the smallest ordinary p of 106 (eighth by HC3 p)' in fig4.md, PAPER and RESULTS.",
    commands=first["R04-04"]["commands"] + "; repair: pandas rank of delta_pvalue and delta_hc3_pvalue by analysis_universe in ambient_program_effects.tsv")
upd["R04-05"] = dict(changed_in_repair=True, panels=P4C_REP,
    what_changes=first["R04-05"]["what_changes"] + " Repair: S4C and 4C re-rendered with forest-axis ticks only inside the plotted range and a full-range axis line; the S4C '-3' label no longer overlaps 'Steatohepatitis'. No value or row changes.",
    commands=first["R04-05"]["commands"] + "; repair: " + REND,
    source_tables=first["R04-05"]["source_tables"] + [f"{FIG}/repair/source_tables/fig4c_h59_registry_117_with_familywise_interval.tsv"],
    note=first["R04-05"]["note"] + " " + INSTALL)
upd["R04-06"] = dict(changed_in_repair=False)
upd["R04-07"] = dict(changed_in_repair=True,
    reproduced=first["R04-07"]["reproduced"] + " Repair: V3 donor_scores_hepatocyte.tsv stages are Healthy 6, MASLD_F0 3, MASH_F2_F3 4, MASH_cirrhosis 4; C2 (disease=1, n=11) therefore includes the 4 cirrhotic donors. H90 (final, untestable) opened GSE289173 and stopped at G3: Ductular injury retains 14 genes and 0.636 of the frozen weight (needs 0.80).",
    what_changes=first["R04-07"]["what_changes"] + " Repair: fig4.md, RESULTS and METHODS state that the C2 disease group includes 4 MASH-cirrhosis donors; fig4.md says 'no test passes Holm correction' instead of 'neither program replicates'; RESULTS calls V3 the only completed new-donor test and adds the H90 supplement sentence (the only text its prespec allows).",
    commands=first["R04-07"]["commands"] + "; repair: pandas value_counts on V3 results/donor_scores_hepatocyte.tsv; read H90 RESULT_H90.md",
    source_tables=first["R04-07"]["source_tables"] + ["/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V3-20261008T015239Z/results/donor_scores_hepatocyte.tsv", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/RESULT_H90.md"],
    note="V3 has no RESULT file. Adding GSE256398 also needs a data-availability and dataset-roster line, not drafted here. The H90 sentence is copied from RESULT_H90.md with its numbers; no H90 number enters a claim.")
upd["R04-08"] = dict(changed_in_repair=False)
upd["R04-09"] = dict(changed_in_repair=False)
upd["R04-10"] = dict(changed_in_repair=True,
    reproduced=first["R04-10"]["reproduced"] + " Repair: stromal_ecm_ductular_by_arm.tsv, family_size 224 in both 774 arms: new774_three_cohorts (364 units) unadjusted fibrosis q 0.214201, composition-adjusted 0.421163. H85 D4_focal_programs.tsv: four-cohort composition-adjusted fibrosis q 0.230276 (rebuilt fractions) and 0.254406 (earlier fractions); H85 PC-3 requires describing ECM only with this range across sets.",
    what_changes="Intro:97-99 no longer calls ECM a 'disease-associated program'. It states that the ECM bulk fibrosis association reaches BH q=0.045 (224 program-by-axis tests) only in the 457-record four-cohort arm without composition adjustment, is q=0.214 in the 364-unit three-cohort arm, and is q=0.23 to 0.25 after composition adjustment depending on the cell-fraction set (H85). 'rests on' is removed.",
    old_value="'a disease-associated program'",
    new_value="bulk fibrosis q 0.045 (457-record arm, unadjusted, 224 tests); 0.214 (364-unit arm); 0.23 to 0.25 composition-adjusted",
    source_tables=first["R04-10"]["source_tables"] + ["/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H85-20261015T124643Z/run-20261015T130425Z/D4_focal_programs.tsv"],
    commands="awk on stromal_ecm_ductular_by_arm.tsv (all arms); column -t on H85 D4_focal_programs.tsv",
    note="The sentence names both arms so it holds whichever arm conflict 2 (not owned here) puts in Figure 4E. Shared line 01_intro:99 with R04-01: apply as string replacements; verified sequentially. The truncated Knapp-Hartung qualification (R12-01) belongs to F4efgh.")
upd["R10-01"] = dict(changed_in_repair=True,
    what_changes=first["R10-01"]["what_changes"] + " Repair: RESULTS:1635 no longer puts 27/113 and 26/114 in one sentence; fig4.md and RESULTS name the arms by record count (457-record four-cohort, 364-unit three-cohort) instead of calling the 457-record arm 'the Figure 4E arm'.",
    new_value="457-record four-cohort arm: 20/114 vs 26.91, P=0.871; NAS 2 vs 2.565, P=0.687. 364-unit three-cohort arm: 1 vs 1.24, P=0.473",
    note=first["R10-01"]["note"] + " Which arm Figure 4E shows is conflict 2 (not owned here); the text is correct under either choice.")
upd["R13-04"] = dict(changed_in_repair=True, blocked_on="",
    what_changes=("Intro:54-56, fig1.md:44-45 and PAPER idea 2 (122-123) stop saying the single-cell programs 'resolve' "
                  "multicellular remodeling; the abstract edit is in R04-01. fig1.md:44-45 now carries the conflict-24 sentence: "
                  "'Figure 4 reports donor-level program scores in five cell lineages and their bulk-tissue context; "
                  "transcription-factor activity is supplementary.' The six conditional title.md rows are withdrawn: clause 2 is "
                  "removed under TITLE_CLAUSE_DECISION.md (row H52 untestable / H54 supports), and the replacement title is drafted "
                  "by Claims R13-01 for title.md, PAPER, README and STATUS."),
    new_value="'score prespecified gene programs in each cell lineage'; fig1: 'program scores in five cell lineages and their bulk-tissue context'; no title edit from this entry",
    note=("Title rows dropped because the H52/H54 outcome is decided and Claims owns the title wording (conflict 0, not owned "
          "here). TF regulons are in S4O (legend 05_figure_legends.md:329-331), which supports 'supplementary'. The bulk "
          "count claim '0/224' from the finding is not reproduced: Fig 4E prints 16 fibrosis and 1 NAS calls."))
upd["G1-01"] = dict(changed_in_repair=True,
    reproduced=first["G1-01"]["reproduced"] + " Repair: 400_mcp/410_prep_atlas_cnmf.py:49-50,130 keeps (sample, cell_type) pairs with at least 30 cells for Hepatocytes, Endothelial cells, Fibroblasts, Macrophages and Cholangiocytes; 'sample' is a GSM library (269 rows in mcp/inputs/donor_metadata.tsv); hotspot_io.py:18 reads those inputs; T cells come from 500_prepare_tcell_subset.py with no minimum (62-roster T-cell minimum 1). fig4d_same_atlas_lineage_specificity.tsv has minimum_cells=50 in all 10 rows.",
    what_changes="METHODS:27-28 now states that the ten localization contrasts use the same 50-cell gate, that the 117-program stage models use the discovery inputs (non-T lineages kept per library only with at least 30 cells; T cells without a minimum), and that a 50-cell gate removes the 62-donor Ductular call (q 0.0376 to 0.0765). Under H59 the call is gone at any gate.",
    new_value="117-program models: 30 cells per library for non-T lineages, no minimum for T cells; ten localization contrasts: 50 cells; 50-cell gate sensitivity q 0.076511, 0/117",
    commands=first["G1-01"]["commands"] + "; repair: sed on 410_prep_atlas_cnmf.py and hotspot_io.py; awk on minimum_cells in fig4d_same_atlas_lineage_specificity.tsv",
    source_tables=first["G1-01"]["source_tables"] + ["/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig4_singlecell_programs/source_tables/current_candidate/fig4d_same_atlas_lineage_specificity.tsv"])
upd["G1-02"] = dict(changed_in_repair=False,
    note=first["G1-02"].get("note", "") + " METHODS:637 holds the end of this old_text and the start of the G1-04 old_text; they do not overlap, so apply both as string replacements (verified sequentially).")
upd["G1-03"] = dict(changed_in_repair=True,
    reproduced=first["G1-03"]["reproduced"] + " Repair: H59 results_primary_family.tsv, smallest T-cell perm_p 0.1947 (module 8), perm_q 0.520384.",
    what_changes=first["G1-03"]["what_changes"] + " Repair: the RESULTS row prints the smallest T-cell permutation p as 0.1947 (q 0.520) instead of 'p >= 0.195'.")
upd["G1-04"] = dict(changed_in_repair=False)

out = []
for fid, e in first.items():
    n = dict(e)
    n.pop("blocked_on", None)
    n.update(upd[fid])
    if not n.get("blocked_on"):
        n.pop("blocked_on", None)
    n["text_edits"] = by_fid.get(fid, [])
    out.append(n)
missing = set(by_fid) - set(first)
assert not missing, missing
(REP / "text_edits/structured_entries_repair.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(len(out), sum(len(o["text_edits"]) for o in out))
for o in out:
    print(o["fid"], o["status"], o["changed_in_repair"], len(o["text_edits"]))
