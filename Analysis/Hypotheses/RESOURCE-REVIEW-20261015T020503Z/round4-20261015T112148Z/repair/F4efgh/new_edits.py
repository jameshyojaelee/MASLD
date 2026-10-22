"""Repair-round text edits for workstream F4efgh (review round 3 -> round 4 repair).

Each edit: (fid, file, line, old_text, new_text, claim_level). old_text is checked
verbatim (exactly one occurrence) against the current file by build_entries.py.
"""

E = []


def add(fid, file, line, old, new, claim=False):
    E.append(dict(fid=fid, file=file, line=str(line), old_text=old, new_text=new, claim_level=claim))


LEG = 'docs/manuscript/05_figure_legends.md'
F4 = 'docs/manuscript/draft/fig4.md'
RES = 'docs/RESULTS.md'
PAP = 'docs/PAPER.md'
STA = 'docs/STATUS.md'
ROAD = 'docs/ROADMAP.md'
MET = 'docs/manuscript/METHODS.md'
F4README = 'figures/main/fig4_singlecell_programs/README.md'
F4CAP = 'figures/main/fig4_singlecell_programs/captions.md'
INDEX = 'figures/main/INDEX.md'

# ---------------------------------------------------------------- R05-02 (4G and 4H on 367 disease units)
add('R05-02', LEG, '268-278',
    "**(G)** Stage-adjusted continuum coefficients overlaid on the fixed,\n"
    "signature-excluded pathway/Hotspot membership geometry and its 43 systems.\n"
    "The 774-unit substrate tests 24,115 genes; memberships, weights and coordinates\n"
    "are unchanged. **(H)** Fibrosis-stage coefficient (score SD per cohort\n"
    "stage-code SD) versus stage-adjusted continuum coefficient (score SD per\n"
    "continuum SD) for all 43 systems, using 435 stage-complete units (109\n"
    "GSE162694, 326 GSE213621). Spearman rho is 0.7265 and 37/43 effects share\n"
    "direction; both summaries are descriptive. One outline identifies the system\n"
    "exceeding the uniform size-matched competitive null on both co-primary axes;\n"
    "none exceeds the connected-region null. The 18 all-four strict maxT calls are\n"
    "a different test and are not outlined in the main panel.",
    "**(G)** Stage-adjusted continuum coefficients overlaid on the fixed,\n"
    "signature-excluded pathway/Hotspot membership geometry and its 43 systems.\n"
    "The 774-unit substrate tests 24,115 genes; memberships, weights and coordinates\n"
    "are unchanged. Coefficients in (G) and (H) are estimated in the 367\n"
    "stage-complete disease units (109 GSE162694, 258 GSE213621); the 68 GSE213621\n"
    "controls have no recorded histology and are excluded from both axes. **(H)**\n"
    "Fibrosis-stage coefficient (score SD per cohort stage-code SD) versus\n"
    "stage-adjusted continuum coefficient (score SD per continuum SD) for all 43\n"
    "systems. Spearman rho is 0.63 and 32/43 effects share direction; both\n"
    "summaries are descriptive. The system maxT and competitive-null tests were run\n"
    "on 435 units that include those controls and are not marked in this panel.")

add('R05-02', LEG, '203',
    "stage-complete anchor models retain 435 genetic-identity clusters.",
    "stage-complete anchor models retain 435 genetic-identity clusters, including\n"
    "the 68 GSE213621 controls, which have no recorded histology and are placed at\n"
    "stage 0; Figures 4G and 4H use the 367 disease units.")

add('R05-02', LEG, '317-319',
    "Matched evaluation-cohort inferential panels compare standardized linear fibrosis-stage\n"
    "trends with stage-adjusted fixed-projection coefficients for the complete\n"
    "transcript, pathway, 117-program Hotspot, and ten-axis NMF families.",
    "Matched evaluation-cohort inferential panels compare standardized linear fibrosis-stage\n"
    "trends with stage-adjusted fixed-projection coefficients for the complete\n"
    "transcript, pathway, 117-program Hotspot, and ten-axis NMF families. These\n"
    "matched comparisons use the 435 stage-complete units, including the 68\n"
    "GSE213621 controls placed at stage 0.")

add('R05-02', F4, '156-160',
    "Stage-adjusted continuum\n"
    "coefficients are overlaid on that geometry. Among all 43 systems, stage and\n"
    "continuum effects have Spearman rho=0.7265 and 37 share direction. The\n"
    "stage-complete evaluation roster contains 435 units, 109 GSE162694 and 326\n"
    "GSE213621.",
    "Stage-adjusted continuum\n"
    "coefficients, estimated in the 367 stage-complete disease units (109 GSE162694\n"
    "and 258 GSE213621), are overlaid on that geometry; the 68 GSE213621 controls\n"
    "have no recorded histology and are excluded. Among all 43 systems, stage and\n"
    "continuum effects have Spearman rho=0.63 and 32 share direction.")

add('R05-02', MET, '371-373',
    "No step producing a\n"
    "Figure 4G or 4H value reads a control label or disease-versus-control result;\n"
    "fixed system memberships and outcome-blind coordinates remain unchanged.",
    "In the five-cohort\n"
    "manifest the GSE213621 stage code places its 68 controls, which have no\n"
    "recorded histology, at stage 0. Figure 4G and 4H stage and continuum\n"
    "coefficients are therefore estimated in disease units only (109 GSE162694 and\n"
    "258 GSE213621). The system maxT and competitive-null tests and the\n"
    "supplementary stage-versus-continuum layer comparisons were run on 435 units\n"
    "in which those controls form their own stage category, and the system\n"
    "redundancy summary on the 466-unit descriptive roster. Fixed system\n"
    "memberships and outcome-blind coordinates remain unchanged.")

add('R05-02', RES, '1362-1363',
    "rho=0.658 and 36/43 share direction (844 records; 0.7265 and 37/43 on the 774\n"
    "units, see Main Figure 4H below)",
    "rho=0.658 and 36/43 share direction (844 records; 0.7265 and 37/43 on the 774\n"
    "units with the 68 GSE213621 controls coded as stage 0; 0.631 and 32/43 on the\n"
    "367 disease units; see Main Figure 4H below)")

add('R05-02', RES, '1268-1269',
    "The stage-adjusted models use 435 units\n"
    "(109 GSE162694 and 326 GSE213621;",
    "The stage-adjusted models use 435 units\n"
    "(109 GSE162694 and 326 GSE213621, the latter including 68 controls that the\n"
    "stage code places at stage 0; the last column of the table below refits both\n"
    "axes without them;")

add('R05-02', RES, '1276-1286',
    "| Layer | Complete family / testable | Continuum-associated | Matched fibrosis-trend comparison |\n"
    "|---|---:|---:|---:|\n"
    "| Ordinary transcripts | 24,115 / 23,972 | 15,215 | 11,130 stage-supported; 7,334 shared; effect rho=0.356 |\n"
    "| Hotspot programs | 117 / 116 | 93 | 64 stage-supported; 54 shared; effect rho=0.823 |\n"
    "| Hallmark | 50 / 45 | 38 | 30 stage-supported; 29 shared; effect rho=0.745 |\n"
    "| KEGG | 658 / 302 | 206 | 145 stage-supported; 126 shared; effect rho=0.669 |\n"
    "| Reactome | 1,787 / 1,124 | 873 | 675 stage-supported; 593 shared; effect rho=0.687 |\n"
    "| GO biological process | 7,583 / 3,317 | 2,567 | 2,163 stage-supported; 1,933 shared; effect rho=0.740 |\n"
    "| GO molecular function | 1,855 / 824 | 539 | 409 stage-supported; 323 shared; effect rho=0.683 |\n"
    "| GO cellular component | 1,042 / 516 | 385 | 277 stage-supported; 242 shared; effect rho=0.700 |\n"
    "| Fixed-basis NMF axes | 10 / 8 with a stage estimate | 8 | 7 stage-supported; 6 shared; effect rho=0.476 |",
    "| Layer | Complete family / testable | Continuum-associated | Matched fibrosis-trend comparison | Effect rho, 367 disease units |\n"
    "|---|---:|---:|---:|---:|\n"
    "| Ordinary transcripts | 24,115 / 23,972 | 15,215 | 11,130 stage-supported; 7,334 shared; effect rho=0.356 | not recomputed |\n"
    "| Hotspot programs | 117 / 116 | 93 | 64 stage-supported; 54 shared; effect rho=0.823 | 0.629 (81/116 same direction) |\n"
    "| Hallmark | 50 / 45 | 38 | 30 stage-supported; 29 shared; effect rho=0.745 | 0.662 |\n"
    "| KEGG | 658 / 302 | 206 | 145 stage-supported; 126 shared; effect rho=0.669 | 0.675 |\n"
    "| Reactome | 1,787 / 1,124 | 873 | 675 stage-supported; 593 shared; effect rho=0.687 | 0.612 |\n"
    "| GO biological process | 7,583 / 3,317 | 2,567 | 2,163 stage-supported; 1,933 shared; effect rho=0.740 | 0.671 |\n"
    "| GO molecular function | 1,855 / 824 | 539 | 409 stage-supported; 323 shared; effect rho=0.683 | 0.523 |\n"
    "| GO cellular component | 1,042 / 516 | 385 | 277 stage-supported; 242 shared; effect rho=0.700 | 0.547 |\n"
    "| Fixed-basis NMF axes | 10 / 8 with a stage estimate | 8 | 7 stage-supported; 6 shared; effect rho=0.476 | not recomputed |\n"
    "\n"
    "The first four columns use the 435 stage-complete units, including the 68\n"
    "GSE213621 controls placed at stage 0. The last column refits both axes in the\n"
    "367 disease units (109 GSE162694, 258 GSE213621; review round 3, R05-02,\n"
    "`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/feature_layer_stage_vs_continuum_rho_by_roster.tsv`).\n"
    "The stage-supported, shared and continuum-associated counts were not\n"
    "recomputed on that roster.")

add('R05-02', RES, '1402',
    "Stage-adjusted n=435 (109 +326); these are",
    "Stage-adjusted n=435 (109 +326, including the 68 GSE213621 controls, which\n"
    "carry no recorded histology); these are")

add('R05-02', RES, '1406',
    "and no step that produces a Figure 4G or 4H value reads a control label or a disease-versus-control result; these values stand on the corrected inputs",
    "but the GSE213621 stage code places its 68 controls at stage 0, so the installed Figure 4G and 4H models include those controls as their own stage category and the 4H stage axis reads the control label (review round 3, R05-02; candidate 4G and 4H on the 367 disease units in `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/repair/fig4g/` and `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`); the other values stand on the corrected inputs")

add('R05-02', RES, '1429-1430',
    "On the 774 units the panel source gives Spearman rho 0.7265\n"
    "(printed as 0.73) and 37/43 in the same direction",
    "On the 774 units the panel source gives Spearman rho 0.7265\n"
    "(printed as 0.73) and 37/43 in the same direction; with the 68 GSE213621\n"
    "controls removed from both axes (367 disease units) the values are 0.631 and\n"
    "32/43 (candidate `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`)")

add('R05-02', PAP, '32-34',
    "Figure 4H carries one outline, on the single system exceeding\n"
    "the uniform size-matched null (0 connected, 1 uniform); the 18 strict maxT\n"
    "systems are not outlined.",
    "Figure 4H plots 43 system coefficients estimated in the 367\n"
    "disease units and carries no outline. On the 435-unit roster that includes the\n"
    "68 GSE213621 controls, one system exceeds the uniform size-matched null and\n"
    "none the connected-region null (0 connected, 1 uniform), and 18 pass strict\n"
    "maxT; none of these is marked in Figure 4H.",
    claim=True)

add('R05-02', PAP, '460-461',
    "Figure 4G places the\n"
    "stage-adjusted continuum coefficients on a signature-excluded, outcome-blind",
    "Figure 4G places the\n"
    "stage-adjusted continuum coefficients, estimated in the 367 disease units, on a\n"
    "signature-excluded, outcome-blind")

add('R05-02', PAP, '474-477',
    "Its single outline marks the one system exceeding 5,000 size-matched\n"
    "draws from the same membership geometry on both co-primary axes. On the 774\n"
    "units, Spearman rho is 0.7265 and 37/43 effects share direction; both remain\n"
    "descriptive.",
    "Both coefficients are estimated in the 367 stage-complete disease units;\n"
    "the 68 GSE213621 controls have no recorded histology. Spearman rho is 0.63\n"
    "and 32/43 effects share direction (0.7265 and 37/43 when the controls entered\n"
    "as stage 0); both remain descriptive. The size-matched competitive null was\n"
    "run on the 435-unit roster with the controls, so the panel carries no outline\n"
    "until it is rerun on disease units.")

add('R05-02', PAP, '484-485',
    "These are\n"
    "different tests; the stage-adjusted roster has 435 units.",
    "These are\n"
    "different tests; their stage-adjusted roster has 435 units, including the 68\n"
    "GSE213621 controls, whereas Figures 4G and 4H use the 367 disease units.")

add('R05-02', PAP, '493-495',
    "Those models\n"
    "   are complete-case in recorded stage and therefore use 109 GSE162694 and 326\n"
    "   GSE213621 units, 435 in total,",
    "Those models\n"
    "   use 109 GSE162694 and 326 GSE213621 units, 435 in total, where the 326\n"
    "   include 68 controls with no recorded histology that the stage code places at\n"
    "   stage 0,")

add('R05-02', STA, '521-522',
    "stage and continuum system effects have Spearman rho\n"
    "0.73 with 37 of 43 the same sign",
    "stage and continuum system effects have Spearman rho\n"
    "0.73 with 37 of 43 the same sign (0.63 and 32 of 43 on the 367 disease units\n"
    "once the 68 GSE213621 controls coded as stage 0 are removed; review round 3,\n"
    "R05-02)")

add('R05-02', STA, '818',
    "Main 4H has one uniform-null outline, rho 0.7265 and 37/43 same direction.",
    "The 326 GSE213621 units include 68 controls coded as stage 0. After review round 3 (R05-02), main 4G and 4H use the 367 disease units; 4H has no outline, rho 0.631 and 32/43 same direction (0.7265 and 37/43 with the controls).")

add('R05-02', STA, '303',
    "main 4H has one uniform-null outline, while the 18 strict maxT calls remain a distinct supplementary criterion.",
    "after review round 3 (R05-02), main 4G and 4H use the 367 disease units and 4H has no outline, while the 18 strict maxT calls (435 units, including the 68 GSE213621 controls) remain a distinct supplementary criterion.")

add('R05-02', ROAD, '282',
    "stage-complete inference 435 units, two co-primary scores; 33 primary BH/direction calls and 18 strict maxT calls distinct from competitive specificity (0 connected, 1 uniform). Main 4H has one uniform-null outline, rho 0.7265 and 37/43 same direction.",
    "stage-complete inference 435 units, including 68 GSE213621 controls with no recorded histology coded as stage 0, two co-primary scores; 33 primary BH/direction calls and 18 strict maxT calls distinct from competitive specificity (0 connected, 1 uniform), all on those 435 units. After review round 3 (R05-02), main 4G and 4H use the 367 disease units; 4H has no outline, rho 0.63 and 32/43 same direction.")

add('R05-02', F4README, '39',
    "| 4G | Stage-adjusted fixed-projection effects across 43 outcome-blind molecular systems | corrected 774-unit panel installed, 2026-10-14; fixed geometry |",
    "| 4G | Stage-adjusted fixed-projection effects across 43 outcome-blind molecular systems | 367 disease units after review round 3 (68 GSE213621 controls excluded); fixed geometry |")

add('R05-02', F4README, '40',
    "| 4H | Fibrosis-stage versus stage-adjusted continuum effects across all 43 systems, with competitive specificity | corrected panel installed, 2026-10-14; one uniform-null outline |",
    "| 4H | Fibrosis-stage versus stage-adjusted continuum effects across all 43 systems | 367 disease units after review round 3; rho 0.63, 32/43 same direction; no outline |")

# ---------------------------------------------------------------- R05-05
add('R05-05', LEG, '299-306',
    "The focal continuum forest shows disease-only stage/sex-adjusted and\n"
    "composition-adjusted effects on the 774-unit disease-only evaluation roster\n"
    "(109 GSE162694 and 258 GSE213621 units). The rebuilt 16-type ILR adjustment\n"
    "reduces the Stromal ECM meta-coefficient from 0.760 to 0.294 (95% CI\n"
    "0.154–0.434) and Ductular injury from 0.707 to 0.601 (0.500–0.701).\n"
    "Stromal ECM is not individually supported after adjustment in GSE162694\n"
    "(0.172, −0.032–0.377; 117-program-family BH q=0.233). Adjustment depends\n"
    "on the deconvolution input; no composition-adjusted system is highlighted.",
    "The focal continuum forest shows stage/sex-adjusted and composition-adjusted\n"
    "effects on the 774-unit disease-only evaluation roster (109 GSE162694 and 258\n"
    "GSE213621 units), with the composition adjustment shown for both cell-type\n"
    "fraction sets. Ductular injury changes little with composition adjustment\n"
    "(meta coefficient 0.60 to 0.73 across codings and fraction sets; unadjusted\n"
    "0.707). Stromal ECM is reduced by composition adjustment by an amount that\n"
    "depends on the deconvolution input (0.203 to 0.479; unadjusted 0.760). In\n"
    "GSE162694 alone, adjusted Stromal ECM is 0.097 to 0.217 with the rebuilt\n"
    "fractions and 0.344 to 0.433 with the earlier fractions. No\n"
    "composition-adjusted program or system is named as a finding.")

add('R05-05', STA, '591-594',
    "composition-adjusted 0.294 [0.154, 0.434], BH q 1.4e-4 (0.726 and 0.488 on the\n"
    "844 records); Ductular injury 0.707 and 0.601 (0.677 and 0.639). Adjusted\n"
    "Stromal ECM in GSE162694 alone is 0.172 [-0.032, 0.377], BH q 0.233, no longer\n"
    "supported there.",
    "composition-adjusted 0.294 [0.154, 0.434], BH q 1.4e-4 in arm A with the rebuilt\n"
    "fractions (0.203 to 0.479 across arms and fraction sets; 0.726 and 0.488 on the\n"
    "844 records); Ductular injury 0.707 and 0.601 (0.60 to 0.73 across arms; 0.677\n"
    "and 0.639 on the 844 records). Adjusted Stromal ECM in GSE162694 alone is 0.172\n"
    "[-0.032, 0.377], BH q 0.233, in arm A with the rebuilt fractions (0.097 to 0.217\n"
    "across arms); with the earlier fractions it is 0.344 to 0.433.")

# ---------------------------------------------------------------- R05-06
add('R05-06', LEG, '312-314',
    "Figure 4H uses one uniform-null outline, as described above. The supplementary\n"
    "inference registry retains 33 two-axis BH/direction calls and 18 all-four\n"
    "strict maxT calls; the jointly powered cross-cohort maxT test gives 33/43.",
    "The supplementary inference registry, fitted on 435 units in which the 68\n"
    "GSE213621 controls form their own stage category, retains 33 two-axis\n"
    "BH/direction calls and 18 strict maxT calls (maxT P<0.05 in both cohorts on\n"
    "both co-primary axes); the jointly powered cross-cohort maxT test gives 33/43.\n"
    "The two co-primary axes give system coefficients correlated at r=0.999, so the\n"
    "strict criterion is in effect a two-cohort conjunction. The 18 strict calls\n"
    "come from models without composition terms. In the 367 disease units, where\n"
    "the unadjusted fit calls all 18 (BH q<0.05 within 43 systems), 4 of the 18 are\n"
    "called in all 11 estimable composition-adjusted fits, 14 lose the call in at\n"
    "least one fit and 9 are called in 5 or fewer.")

add('R05-06', PAP, '511-512',
    "S4U retains\n"
    "   all 43 systems and marks only those 18 strict calls with a neutral asterisk;",
    "S4U retains\n"
    "   all 43 systems and marks only those 18 strict calls with a neutral asterisk\n"
    "   (in the 367 disease units, 4 of the 18 are called in all 11 estimable\n"
    "   composition-adjusted fits, 14 lose the call in at least one, and 9 are\n"
    "   called in 5 or fewer);")

add('R05-06', F4, '186-193',
    "Figure 4H compares the stage coefficient in score SD per cohort SD of the\n"
    "recorded stage code with the stage/sex-adjusted continuum coefficient in score\n"
    "SD per continuum SD. These estimate different quantities. Its single outline\n"
    "marks the system exceeding the uniform size-matched competitive null on both\n"
    "co-primary axes; no system exceeds the connected-region null on both axes.\n"
    "The separate 18 all-four strict maxT calls remain in the supplement and do not\n"
    "establish competitive specificity. The jointly powered cross-cohort maxT\n"
    "test gives 33 of 43 systems on the corrected counts.",
    "Figure 4H compares the fibrosis-stage coefficient in score SD per cohort SD\n"
    "of the recorded stage code with the stage/sex-adjusted continuum coefficient in\n"
    "score SD per continuum SD, both estimated in the 367 disease units. These\n"
    "estimate different quantities. The system maxT and competitive-null tests were\n"
    "run on 435 units that included the 68 GSE213621 controls. On that roster no\n"
    "system exceeds the connected-region null on both axes, one (S01) exceeds the\n"
    "permissive uniform size-matched null, and 18 pass maxT in both cohorts on both\n"
    "co-primary axes. The two axes give system coefficients correlated at r=0.999,\n"
    "so this is in effect a two-cohort criterion. In the 367 disease units, 14 of\n"
    "the 18 lose the call in at least one of 11 estimable composition-adjusted fits\n"
    "and 9 are called in 5 or fewer. The jointly powered cross-cohort maxT test\n"
    "gives 33 of 43. Because the system scores are strongly correlated and none\n"
    "exceeds the connected-region null, these counts are reported as part of a\n"
    "diffuse shared axis, not as 33 or 18 separate findings.")

# ---------------------------------------------------------------- R12-01
add('R12-01', RES, '1713-1715',
    "The loss of calls in the three-cohort arm follows from the cohort count: the\n"
    "old records without GSE174478 (376 records, three cohorts) give 1 fibrosis and\n"
    "0 NAS calls.",
    "The four-cohort calls are not explained by the cohort count alone. All 20\n"
    "(16 fibrosis, 1 NAS, 3 composition-adjusted fibrosis with the rebuilt\n"
    "fractions) have tau2 = 0 and between-cohort Q of 0.09 to 1.61, below the\n"
    "expected value of 3, so their untruncated Knapp-Hartung SE is 0.17 to 0.73 of\n"
    "the REML SE. With the Knapp-Hartung variance floored at the REML value\n"
    "(metafor test=\"adhoc\"), both arms give 0 calls in every model (smallest q\n"
    "0.193 in the four-cohort arm). A REML z-test gives 47/28/21/10 calls in the\n"
    "four-cohort arm and 37/14/13/11 in the three-cohort arm; with three or four\n"
    "cohorts it ignores the uncertainty in tau2, and the Figure 4E code uses\n"
    "Knapp-Hartung (`scripts/manuscript/program_observability_map/analysis_lib.R`)\n"
    "(`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4e/source_tables/fig4e_calls_by_arm_and_hk_variant.tsv`).\n"
    "With the earlier cell-type fractions the four-cohort composition-adjusted\n"
    "count is 1 fibrosis and 0 NAS instead of 3 and 0 (H85), and the three-cohort\n"
    "arm stays at 0 and 0 under both fraction sets. The old records without\n"
    "GSE174478 (376 records, three cohorts) give 1 fibrosis and 0 NAS calls.")

# ---------------------------------------------------------------- R13-03 and conflict 2 (Figure 4E arm everywhere)
add('R13-03', LEG, '248-258',
    "**(E)** Weighted program-score linear associations with recorded fibrosis\n"
    "category and NAS in the four-cohort transport arm: 457 records comprising\n"
    "364 selected units from GSE130970, GSE135251 and GSE162694 plus 93 source\n"
    "records from GSE174478. This arm is separate from the five-cohort 774-unit\n"
    "disease contrast. Effects are score SD per stage-code unit or NAS point;\n"
    "color shows BH q over 224 tests (112 eligible programs × two axes).\n"
    "Sixteen fibrosis and one NAS association meet BH q <0.05; after composition\n"
    "adjustment, three fibrosis and zero NAS associations remain. GSE174478's\n"
    "source-record identity qualification remains separate from the selected\n"
    "genetic clusters. Linear category coding does not establish nonlinear or\n"
    "longitudinal remodeling.",
    "**(E)** Weighted program-score linear associations with recorded fibrosis\n"
    "category and NAS on the 774-unit substrate: 364 selected genetic-identity\n"
    "clusters from GSE130970, GSE135251 and GSE162694, the three Resource cohorts\n"
    "that record both axes. Each point is one of 112 eligible programs. The x axis\n"
    "shows the effect in score SD per stage-code unit or NAS point from a\n"
    "three-cohort random-effects meta-analysis (REML, Knapp-Hartung); the y axis\n"
    "shows -log10 BH q over 224 tests (112 programs × two axes), and the dotted\n"
    "line marks q = 0.05. No fibrosis or NAS association meets BH q <0.05\n"
    "(smallest q 0.079 and 0.159), with or without composition adjustment under\n"
    "either cell-type fraction set. Adding 93 GSE174478 source records, which lie\n"
    "outside the 774-unit substrate and are not identity-resolved, gives 16\n"
    "fibrosis and one NAS association at q <0.05; after composition adjustment,\n"
    "three fibrosis associations remain with the rebuilt fractions and one with the\n"
    "earlier fractions. Each of these calls has between-cohort Q below its expected\n"
    "value, and none remains when the Knapp-Hartung variance is not allowed to fall\n"
    "below the REML random-effects variance. Linear category coding does not\n"
    "establish nonlinear or longitudinal remodeling.")

add('R13-03', F4, '5-6',
    "Figure 4E separately\n"
    "> includes 93 source records from GSE174478 outside that substrate.",
    "Figure 4E uses the\n"
    "> 774-unit substrate only; 93 GSE174478 source records enter a sensitivity arm.")

add('R13-03', F4, '113-121',
    "Figure 4E shows linear fibrosis-category and NAS associations for all 112\n"
    "eligible weighted programs, with BH correction over 224 tests. The selected\n"
    "four-cohort arm uses 364 units from three corrected bulk cohorts plus 93\n"
    "GSE174478 source records, 457 in total. Sixteen fibrosis and one NAS association\n"
    "meet BH q<0.05; after composition adjustment, three fibrosis and no NAS\n"
    "association remain. The arm differs from the five-cohort pooled disease\n"
    "contrast and retains GSE174478's source-record identity limitation. Linear\n"
    "category coding does not summarize nonlinear remodeling or within-person\n"
    "change.",
    "Figure 4E shows linear fibrosis-category and NAS associations for all 112\n"
    "eligible weighted programs, with BH correction over 224 tests, in the 364\n"
    "selected genetic-identity clusters of the three 774-unit cohorts that record\n"
    "both fibrosis and NAS. No program meets BH q<0.05 on either axis (smallest q\n"
    "0.079 for fibrosis and 0.159 for NAS), with or without composition adjustment\n"
    "under either cell-type fraction set. Adding 93 GSE174478 source records, which\n"
    "are outside the 774-unit substrate and not identity-resolved, raises the\n"
    "meta-analysis from three to four cohorts and gives 16 fibrosis and one NAS\n"
    "call; after composition adjustment three fibrosis calls remain with the\n"
    "rebuilt fractions and one with the earlier fractions. All of these calls\n"
    "depend on the untruncated Knapp-Hartung variance: each has between-cohort Q\n"
    "below its expected value, and none remains when that variance is floored at\n"
    "the REML value. Figure 4E therefore shows no program association with linear\n"
    "fibrosis category or NAS that holds across these analysis choices. Linear\n"
    "category coding does not summarize nonlinear remodeling or within-person\n"
    "change.")

add('R13-03', PAP, '444-450',
    "The selected four-cohort arm uses GSE130970, GSE135251 and GSE162694\n"
    "(364 selected genetic-identity units) plus GSE174478 (93 source records),\n"
    "457 in total. GSE174478 is outside the 774-unit pooled disease fit and\n"
    "retains a separate source-record identity limitation. Sixteen fibrosis and\n"
    "one NAS association meet BH q<0.05; after composition adjustment, three\n"
    "fibrosis and zero NAS associations remain. The three-cohort sensitivity\n"
    "has no call.",
    "Figure 4E uses GSE130970, GSE135251 and GSE162694 (364 selected\n"
    "genetic-identity units), the three 774-unit cohorts that record fibrosis and\n"
    "NAS. No fibrosis or NAS association meets BH q<0.05 (smallest q 0.079 and\n"
    "0.159), with or without composition adjustment under either fraction set. A\n"
    "sensitivity arm that adds 93 GSE174478 source records (outside the 774-unit\n"
    "fit, not identity-resolved) gives 16 fibrosis and one NAS call; after\n"
    "composition adjustment three fibrosis calls remain with the rebuilt fractions\n"
    "and one with the earlier fractions. All depend on the untruncated\n"
    "Knapp-Hartung variance (Q below its expectation; none under the truncated\n"
    "variant) and are not reported as findings.")

add('R13-03', RES, '1727-1729',
    "The 2026-10-13 Figure 4 candidate uses the four-cohort arm (774 units plus\n"
    "GSE174478, named in the legend), which the main session recommended and the\n"
    "user told it to proceed with. The three-cohort version is kept beside it.",
    "The 2026-10-13 Figure 4 candidate used the four-cohort arm (774 units plus\n"
    "GSE174478, named in the legend), which the main session recommended and the\n"
    "user told it to proceed with; both arms had been computed on 2026-10-12 before\n"
    "that choice. After review round 3 (R13-03, R12-01) Figure 4E uses the\n"
    "three-cohort 774-unit arm (0 fibrosis and 0 NAS calls), and the four-cohort\n"
    "arm is a sensitivity whose calls need GSE174478 and the untruncated\n"
    "Knapp-Hartung variance (candidate\n"
    "`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4e/`).")

add('R13-03', STA, '608-613',
    "On 2026-10-13 the user told the main session to proceed with its\n"
    "recommendations: Figure 4E is the four-cohort arm (774 units plus GSE174478,\n"
    "named in the legend), and Figure 4H carries a single outline on the one system\n"
    "exceeding the size-matched null, with no maxT count printed (`docs/PAPER.md`\n"
    "now says this in both passages). The other 4E and 4H versions are kept in the\n"
    "same directory.",
    "On 2026-10-13 the user told the main session to proceed with its\n"
    "recommendations: Figure 4E is the four-cohort arm (774 units plus GSE174478,\n"
    "named in the legend), and Figure 4H carries a single outline on the one system\n"
    "exceeding the size-matched null, with no maxT count printed (`docs/PAPER.md`\n"
    "now says this in both passages). The other 4E and 4H versions are kept in the\n"
    "same directory. Review round 3 replaced both choices. Figure 4E uses the\n"
    "three-cohort 774-unit arm (0 and 0 calls), because all 20 four-cohort calls\n"
    "need GSE174478 and the untruncated Knapp-Hartung variance (R13-03, R12-01).\n"
    "Figures 4G and 4H use the 367 disease units, with no outline in 4H, because\n"
    "the GSE213621 stage code placed its 68 controls at stage 0 (R05-02).")

add('R13-03', STA, '6-7',
    "Figure 4A is complete, 4E uses the selected 457-record transport arm, and 4H\n"
    "has one competitive-null outline.",
    "Figure 4A is complete; after review round 3, 4E uses the three-cohort 364-unit\n"
    "arm (the 457-record arm is a sensitivity), and 4G and 4H use the 367 disease\n"
    "units, with no outline in 4H.")

add('R13-03', STA, '585-587',
    "composition adjustment; 774 units only (three cohorts, 364 units, 112 eligible)\n"
    "0 and 0 calls; 774 units plus GSE174478 (457 records, 112 eligible) 16 and 1\n"
    "calls, 3 and 0 after composition adjustment.",
    "composition adjustment; 774 units only (three cohorts, 364 units, 112 eligible;\n"
    "Figure 4E after review round 3) 0 and 0 calls, also after composition\n"
    "adjustment with either fraction set; 774 units plus GSE174478 (457 records, 112\n"
    "eligible; sensitivity) 16 and 1 calls, 3 and 0 after composition adjustment\n"
    "with the rebuilt fractions and 1 and 0 with the earlier fractions (H85), and 0\n"
    "and 0 under the truncated Knapp-Hartung variance (R12-01).")

add('R13-03', STA, '814',
    "selected Figure 4E arm is four-cohort 457 records, separate from the pooled 774-unit fit",
    "after review round 3, Figure 4E uses the three-cohort 364-unit arm of the 774-unit substrate; the four-cohort 457-record arm, which adds GSE174478 from outside the pooled 774-unit fit, is a sensitivity")

add('R13-03', STA, '888',
    "4E uses 112 programs in 457 records, 4H one uniform-null outline;",
    "4E uses 112 programs in the three-cohort 364-unit arm (the 457-record arm is a sensitivity), 4G and 4H use the 367 disease units with no outline in 4H (review round 3);")

add('R13-03', 'docs/README.md', '6-7',
    "Figure 4 includes the atlas schematic,\n"
    "457-record transport arm and one competitive-null outline.",
    "Figure 4 includes the atlas schematic,\n"
    "the three-cohort 364-unit Figure 4E arm (review round 3) and Figure 4G/H on\n"
    "the 367 disease units, with no outline in 4H.")

add('R13-03', ROAD, '299',
    "E uses the selected 457-record transport arm.",
    "E uses the three-cohort 364-unit arm of the 774-unit substrate after review\n"
    "round 3; the 457-record arm is a sensitivity.")

add('R13-03', F4README, '8-12',
    "Figure 4E uses the four-cohort 457-record arm (112 programs, 224-test BH family;\n"
    "16 fibrosis and one NAS call), including 93 GSE174478 source records outside\n"
    "the corrected pooled 774-unit fit. Figure 4H prints rho=0.73 and 37/43 same\n"
    "direction with one uniform size-matched-null outline, not the 18 strict maxT\n"
    "outlines.",
    "After review round 3, Figure 4E uses the three-cohort 364-unit arm of the\n"
    "774-unit substrate (112 programs, 224-test BH family; no fibrosis or NAS call\n"
    "at q<0.05); the four-cohort 457-record arm, which adds 93 GSE174478 source\n"
    "records (16 fibrosis and one NAS call), is a sensitivity. Figures 4G and 4H use\n"
    "the 367 stage-complete disease units; Figure 4H prints rho=0.63 and 32/43 same\n"
    "direction with no outline.")

add('R13-03', F4README, '37',
    "| 4E | Complete eligible weighted linear bulk associations with fibrosis and NAS | corrected four-cohort 457-record arm installed, 2026-10-14 |",
    "| 4E | Complete eligible weighted linear bulk associations with fibrosis and NAS | three-cohort 364-unit arm after review round 3; four-cohort 457-record arm is a sensitivity |")

add('R13-03', F4CAP, '8-11',
    "Figure 4E uses 457 records and 112 eligible programs in a 224-test BH\n"
    "family. Figure 4H shows one uniform size-matched-null outline, descriptive\n"
    "rho=0.7265 and 37/43 effects with the same sign; its 18 strict maxT calls are\n"
    "not outlined.",
    "Figure 4E uses the 364 units of the three 774-unit cohorts that record\n"
    "fibrosis and NAS and 112 eligible programs in a 224-test BH family; no program\n"
    "meets q<0.05. Figures 4G and 4H use the 367 stage-complete disease units;\n"
    "Figure 4H shows descriptive rho=0.63 and 32/43 effects with the same sign and\n"
    "no outline.")

add('R13-03', INDEX, '38-43',
    "Figure 4E uses\n"
    "112 programs in the selected four-cohort 457-record arm, with 16 fibrosis and\n"
    "one NAS call (BH over 224 tests); GSE174478 is outside the pooled 774-unit\n"
    "fit. Figure 4H prints rho=0.73 and 37/43 same direction, with one outline\n"
    "for the uniform size-matched competitive null. Eighteen strict maxT calls\n"
    "remain a separate supplementary criterion.",
    "Figure 4E uses\n"
    "112 programs in the three-cohort 364-unit arm of the 774-unit substrate, with\n"
    "no fibrosis or NAS call (BH over 224 tests); the four-cohort 457-record arm,\n"
    "which adds GSE174478 from outside the pooled 774-unit fit, is a sensitivity.\n"
    "Figures 4G and 4H use the 367 stage-complete disease units; Figure 4H prints\n"
    "rho=0.63 and 32/43 same direction, with no outline. Eighteen strict maxT calls\n"
    "(435 units, including the 68 GSE213621 controls) remain a separate\n"
    "supplementary criterion.")

add('R13-03', F4, '193-195',
    "Likewise, the existing Hotspot matched-background comparison\n"
    "does not establish exceptional disease specificity (fibrosis 26/114 versus\n"
    "null mean 18.445, nominal P=0.134; NAS 4 versus 4.225, P=0.512).",
    "Likewise, the Hotspot matched-background comparison in the\n"
    "364-unit Figure 4E arm does not establish disease specificity: with unweighted\n"
    "scores, 1 of 114 programs reaches fibrosis BH q<0.05 (228-test family)\n"
    "against a mean of 1.24 for size- and expression-matched random gene sets\n"
    "(empirical P=0.473, 200 draws), and none reaches NAS against 0.19 (P=1). In\n"
    "the 457-record sensitivity arm the counts are 20 against 26.9 (P=0.871) and 2\n"
    "against 2.57 (P=0.687).")

add('R13-03', RES, '1635-1636',
    "The 27/113 belongs to the 469-record working roster;\n"
    "the 774-unit arms are in the Figure 4E paragraph below.",
    "The 27/113 and 26/114 belong to the 469-record working\n"
    "roster. On the 364-unit three-cohort Figure 4E arm the unweighted vocabulary\n"
    "count is 1/114 against a matched-random mean of 1.24 (p = 0.473), and NAS 0\n"
    "against 0.19 (p = 1); on the 457-record four-cohort sensitivity arm it is\n"
    "20/114 against 26.91 (p = 0.871), and NAS 2 against 2.565 (p = 0.687). The\n"
    "774-unit weighted arms are in the Figure 4E paragraph below.")

add('R13-03', RES, '2146-2147',
    "the vocabulary null (26 of 114\n"
    "fibrosis-supported against 18.4 expected, p = 0.134) now has a cause,",
    "the vocabulary null (469-record roster: 26 of 114\n"
    "fibrosis-supported against 18.4 expected, p = 0.134; 364-unit Figure 4E arm:\n"
    "1 against 1.24, p = 0.473; 457-record sensitivity arm: 20 against 26.9,\n"
    "p = 0.871) now has a cause,")

add('R13-03', 'docs/manuscript/draft/01_intro.md', '97-99',
    "ECM/IGFBP7 and ductular-injury/BICC1 test whether a disease-associated program,\n"
    "its cell-of-origin assignment, and its intact-tissue organization support the\n"
    "same conclusion.",
    "ECM/IGFBP7 and ductular-injury/BICC1 test whether a program's disease\n"
    "association, its cell-of-origin assignment, and its intact-tissue organization\n"
    "support the same conclusion. For ECM/IGFBP7 the bulk fibrosis association does\n"
    "not reach BH q<0.05 in the 364 units of the three 774-unit cohorts (q=0.214;\n"
    "q=0.421 after composition adjustment); it reaches q=0.045 only when 93\n"
    "GSE174478 source records are added.")

# ---------------------------------------------------------------- R13-05
add('R13-05', 'docs/manuscript/draft/abstract.md', '42-46',
    "An externally anchored molecular continuum\n"
    "complemented recorded fibrosis stage, preserving the major metabolic,\n"
    "inflammatory, ductular, and extracellular-matrix patterns while resolving\n"
    "coordinated transcriptional variation among participants assigned the same\n"
    "stage.",
    "An externally anchored molecular continuum\n"
    "complemented recorded fibrosis stage, preserving the major metabolic,\n"
    "inflammatory, ductular, and extracellular-matrix patterns. Among samples\n"
    "sharing a recorded stage, many program and pathway scores still varied with\n"
    "the continuum; which programs carry that variation depends on adjustment for\n"
    "cell-type composition.",
    claim=True)

add('R13-05', 'docs/manuscript/draft/fig3.md', '88-91',
    "The\n"
    "continuum recovered the broad stage-associated direction of transcript,\n"
    "pathway, NMF, and Hotspot remodeling while retaining substantial stage-adjusted\n"
    "signal, particularly for multicellular programs.",
    "The\n"
    "continuum recovered the broad stage-associated direction of transcript,\n"
    "pathway, NMF, and Hotspot remodeling while retaining stage-adjusted signal.\n"
    "For the Hotspot programs that signal depends on cell-type composition\n"
    "adjustment: in the 367 disease units, 94 of 116 estimable programs are\n"
    "continuum-associated (BH q<0.05) without adjustment and 53 to 76 after\n"
    "adjustment, depending on the adjustment arm and the deconvolution input.")

add('R13-05', F4, '126-128',
    "In GSE162694 and GSE213621, both Stromal ECM and Ductular injury increase\n"
    "along the fixed external projection in stage- and sex-adjusted donor-level\n"
    "models, with concordant effects for the cohort-PC1 score.",
    "In the 367 disease units of GSE162694 and GSE213621, both Stromal ECM and\n"
    "Ductular injury increase along the fixed external projection in stage- and\n"
    "sex-adjusted genetic-cluster-level models (meta coefficients 0.760 and 0.707\n"
    "score SD per continuum SD), with concordant effects for the cohort-PC1 score.\n"
    "After cell-type composition adjustment, Ductular injury changes little (0.60\n"
    "to 0.73 across adjustments and fraction sets), whereas Stromal ECM falls to\n"
    "0.203 to 0.479 depending on the deconvolution input; in GSE162694 alone it is\n"
    "not supported with the rebuilt fractions (0.172 [-0.032, 0.377], BH q 0.233).")

add('R13-05', F4, '139-141',
    "histology\n"
    "anchors severity, whereas the continuum resolves coordinated molecular\n"
    "heterogeneity among participants sharing a recorded histologic category.",
    "histology\n"
    "anchors severity, whereas the continuum describes a shared axis of molecular\n"
    "variation among samples sharing a recorded histologic category.")

add('R13-05', F4, '149-150',
    "auditable groupings so that remodeling within a recorded stage can be described\n"
    "as coordinated rather than gene by gene.",
    "auditable groupings so that remodeling within a recorded stage can be described\n"
    "by feature group rather than gene by gene. The groupings are not independent:\n"
    "the 43 system scores behave like about 14 independent systems, and their first\n"
    "principal component carries 41% of the variance (466-unit descriptive\n"
    "roster).")

add('R13-05', 'docs/manuscript/draft/fig5_discussion.md', '200-203',
    "Its strong alignment with fibrosis-based pathway, NMF, and Hotspot\n"
    "results shows that the major biological programs are preserved, while\n"
    "stage-adjusted associations show that histologic bins compress coordinated\n"
    "molecular variation.",
    "Its strong alignment with fibrosis-based pathway, NMF, and Hotspot\n"
    "results shows that the major biological programs are preserved, while\n"
    "stage-adjusted associations show that molecular variation within histologic\n"
    "bins tracks the continuum; which programs carry that variation depends on\n"
    "adjustment for cell-type composition.")

add('R13-05', 'docs/manuscript/draft/fig5_discussion.md', '207-208',
    "This argues for coordinated yet cohort-variable remodeling rather than\n"
    "a universal molecular sequence or ontology.",
    "This argues for shared yet cohort-variable remodeling rather than\n"
    "a universal molecular sequence or ontology. Program-pair covariance seen in\n"
    "the discovery cohorts recurred in two further bulk cohorts (GSE281797 and\n"
    "GSE225740); their participants were not checked by genotype against the\n"
    "discovery cohorts, and bulk covariance mixes cell abundance with program\n"
    "activity.",
    claim=True)

# ---------------------------------------------------------------- R13-11
add('R13-11', 'docs/manuscript/draft/fig1.md', '44-45',
    "Figure 4\n"
    "resolves donor-level multicellular programs and regulatory activity.",
    "Figure 4\n"
    "resolves donor-level single-cell programs and their bulk-tissue context;\n"
    "transcription-factor activity is supplementary.")

add('R13-11', F4, '1',
    "## Donor-level single-cell programs resolve cell-state and regulatory context",
    "## Donor-level single-cell programs and their bulk-tissue context")

# ---------------------------------------------------------------- unchanged first-round edits carried by fid
UNCHANGED_FROM_ROUND1 = {
    'R05-07': None,  # RESULTS 1416
    'R05-09': None,  # RESULTS 1424
}
