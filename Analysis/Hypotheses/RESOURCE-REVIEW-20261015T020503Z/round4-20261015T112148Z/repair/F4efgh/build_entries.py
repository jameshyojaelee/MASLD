"""Assemble the F4efgh repair output: all 11 entries, checker resolutions and owned conflicts [2, 13, 14].

Starts from the round-1 entries in repair_inputs/F4efgh.json, replaces text_edits with the
validated set (new_edits.py plus the unchanged R05-07 and R05-09 edits), and rewrites the
fields that changed in repair. Runs check_edits.py logic first and refuses to write on any failure.
"""
import copy, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
R4 = ROOT + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/'
CAND = ROOT + 'figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/'
REP = CAND + 'repair/fig4g/'

chk = subprocess.run([sys.executable, HERE + '/check_edits.py'], capture_output=True, text=True,
                     env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
print(chk.stdout)
assert chk.returncode == 0 and 'bad 0' in chk.stdout, chk.stderr
edits = json.load(open(HERE + '/tmp/all_edits.json'))

r1 = {e['fid']: e for e in json.load(open(R4 + 'repair_inputs/F4efgh.json'))['entries']}
out = {}
for fid, e in r1.items():
    e = copy.deepcopy(e)
    e['text_edits'] = [{k: t[k] for k in ('file', 'line', 'old_text', 'new_text', 'claim_level')}
                       for t in edits if t['fid'] == fid]
    if not e['text_edits']:
        e.pop('text_edits')
    e['changed_in_repair'] = False
    out[fid] = e

G4_PDF = REP + 'fig4g_molecular_systems_continuum_effect_map_disease_only.pdf'
REPAIR_CMDS = ('Repair, login node (light render on saved tables, no fitting): rnaseq python '
               'repair/F4efgh/build_fig4g_disease_only_nodes.py (merges the installed 4G node overlay with '
               'round-1 job 7070 betas by roster); rnaseq Rscript repair/F4efgh/render_fig4g_from_nodes.R '
               '(copy of render_G_774.R, 220x176 pt, seed 20260817). The as_published node table re-renders '
               'to a PDF whose pdftotext is identical to the installed 4G (repair/F4efgh/fig4g_reproduction/). '
               'pdftotext and pdftoppm checks of the disease-only PDF. check_edits.py: every old_text found '
               'once in the current file, no overlap among F4efgh edits, no em dash or evidence tag.')

# ---------------------------------------------------------------- R05-02
e = out['R05-02']
e['changed_in_repair'] = True
e['reproduced'] = e['reproduced'] + (
    ' Repair: Figure 4G re-rendered on the same 367 units. The as_published node table reproduces the '
    'installed node betas to 1.3e-15 and its PDF text is identical to the installed panel; on 367 units the '
    'system continuum betas equal the 4H candidate source to 8e-17, S01 is 0.665 (0.712 on 435), and the '
    'colour limit (98th percentile of |beta|) is 0.796 (0.838).')
e['what_changes'] = (
    "Fig 4H printed 'Spearman rho = 0.73; 37/43 same direction' becomes 0.63 and 32/43 on 367 disease units, "
    "with no outline; the x axis becomes a true fibrosis-stage coefficient. Fig 4G is re-rendered on the same "
    "367 units, so one system no longer has two continuum coefficients in one figure (S01 0.665 in both 4G and 4H; "
    "it was 0.712 in the installed 4G); the colour scale runs to 0.796 instead of 0.838 and the geometry is unchanged. "
    "The legend states the 367-unit roster once for (G) and (H) and says the maxT and competitive-null tests are on "
    "435 units with the controls. METHODS:371 ('no 4G/4H value reads a control label') is false and is rewritten. "
    "The S3 stage-versus-continuum layer table gains a 367-unit effect rho column (Hotspot 0.823 to 0.629, 81/116; "
    "Hallmark 0.745 to 0.662; KEGG 0.669 to 0.675; Reactome 0.687 to 0.612; GO BP 0.740 to 0.671; GO CC 0.700 to "
    "0.547; GO MF 0.683 to 0.523). The S3 (I-L) and matched-panel legends now say their 435 units include the 68 "
    "controls at stage 0. S4J is not re-rendered.")
e['panels'] = e['panels'] + [G4_PDF]
e['source_tables'] = e['source_tables'] + [REP + 'fig4g_nodes_disease_only.tsv.gz',
                                           REP + 'fig4g_nodes_as_published.tsv.gz',
                                           REP + 'fig4g_node_checks.tsv']
e['commands'] = e['commands'] + ' ' + REPAIR_CMDS
e['new_value'] = ('4H rho 0.63 (0.631078), 32/43, n=367 (109+258), no outline; 4G on the same 367 units '
                  '(S01 system continuum beta 0.665; colour limit 0.796); S3 layer effect rho on 367: Hotspot '
                  '0.629 (81/116), Hallmark 0.662, KEGG 0.675, Reactome 0.612, GO BP 0.671, GO CC 0.547, GO MF 0.523')
e['old_value'] = ('4H rho 0.73 (0.726518), 37/43, n=435 (109+326), one uniform-null outline; 4G on 435 '
                  '(S01 0.712; colour limit 0.838); S3 layer effect rho Hotspot 0.823 (93/116), Hallmark 0.745, '
                  'KEGG 0.669, Reactome 0.687, GO BP 0.740, GO CC 0.700, GO MF 0.683')
e['note'] = (
    'Not run: the system maxT registry (18 strict, 33 joint), the competitive nulls (1 uniform, 0 connected) and '
    'the effective-system count (14) remain on the 435-unit roster (14 on the 466-unit descriptive roster); '
    '70_infer_molecular_systems.R and 73_stress_test_molecular_systems.R would need a disease-only rerun before '
    'any outline or count returns to 4H. The stage-supported, shared and continuum-associated counts in the S3 '
    'layer table and the transcript and NMF effect rho were not recomputed on 367 units; no running plan covers '
    'them (H51 tests a different estimand). S4J is not re-rendered. 4G disease-only nodes keep the installed '
    'geometry; q and support columns are blanked because they belong to the 435-unit fits and are not drawn. '
    'The combined (G)+(H) legend edit replaces the round-1 (H)-only edit; the fig4.md:156 edit now covers the '
    '4G sentence. Install 4G and 4H candidates together with these legend lines. Conflict 13 and 14 edits '
    '(ROADMAP 282, fig4 README 39-40, STATUS 303, legend 203 and 317, RESULTS layer table) are listed here.')

# ---------------------------------------------------------------- R05-04
e = out['R05-04']
e['changed_in_repair'] = True
e['new_value'] = ('GSE213621 n=258; GSE162694 adjusted Stromal ECM 0.172 [-0.032, 0.377] (rebuilt fractions), '
                  '0.433 [0.228, 0.638] (earlier fractions)')
e['note'] = (e.get('note', '') + ' Repair: the earlier-fraction upper bound is 0.63850 in '
             'meta_effects_all_arms.tsv and prints as 0.638 (was 0.639). No panel or text edit uses it.').strip()

# ---------------------------------------------------------------- R05-05
e = out['R05-05']
e['changed_in_repair'] = True
e['what_changes'] = (
    "The S3-S4 legend stops reporting single-arm adjusted values as the result. It states that Ductular injury "
    "changes little with composition adjustment (0.60 to 0.73; unadjusted 0.707), that Stromal ECM falls to 0.203 "
    "to 0.479 depending on the deconvolution input (unadjusted 0.760), and that GSE162694 Stromal ECM is 0.097 to "
    "0.217 (rebuilt) versus 0.344 to 0.433 (earlier). 'Holds' is replaced by 'changes little' because no adjusted "
    "program is named as a finding (H85 reading rule; H59 removes Ductular injury as a primary single-cell call). "
    "STATUS (5) gets the same ranges.")
e['note'] = ("RESULTS:285 records the user's 2026-10-12 reporting rule, which uses 'holding'; it is a record of "
             "that decision and is left unchanged.")

# ---------------------------------------------------------------- R05-06
e = out['R05-06']
e['changed_in_repair'] = True
e['what_changes'] = (
    "No number changes. The S3-S4 legend, the fig4.md 4H paragraph and the PAPER S4U sentence say that the 18 "
    "strict calls come from the 435-unit roster with controls and from models without composition terms, and that "
    "in the 367 disease units (where the unadjusted fit calls all 18) 4 of the 18 are called in all 11 estimable "
    "composition-adjusted fits, 14 lose the call in at least one and 9 are called in 5 or fewer. The overstatement "
    "'most S4U asterisks mark systems whose call does not survive adjustment' and the '14 effective systems, one "
    "shared axis' inference are removed; the shared-axis reading now rests on the strong correlation of system "
    "scores and 0/43 connected-region null. S4U is not re-rendered.")
e['new_value'] = ('18 strict maxT calls (435-unit roster, no composition terms); on 367 units 4/18 called in all 11 '
                  'estimable composition-adjusted fits, 14/18 lose the call in at least one, 9/18 called in 5 or '
                  'fewer (distribution 0:1, 2:2, 3:2, 4:1, 5:3, 6:1, 10:4, 11:4)')
e['commands'] = e['commands'] + ' Repair: per-system call counts over the 11 fits tabulated from the same join.'

# ---------------------------------------------------------------- R12-01
e = out['R12-01']
e['changed_in_repair'] = True
e['what_changes'] = e['what_changes'] + (
    " Repair: the RESULTS text says 'both arms' (the four-cohort arm is outside the 774 units), names "
    "Knapp-Hartung as the method coded for Figure 4E instead of calling the z-test 'not prespecified', and adds "
    "the H85 earlier-fraction count (1 fibrosis, 0 NAS after composition in the four-cohort arm; the three-cohort "
    "arm stays 0 and 0 under both fraction sets).")
e['new_value'] = e['new_value'] + '; four-cohort composition-adjusted 1/0 with the earlier fractions (H85), 3/0 rebuilt'
e['source_tables'] = e['source_tables'] + [ROOT + 'Analysis/Hypotheses/H85-20261015T124643Z/RESULT_H85.md']
e['note'] = (e['note'] + " knha is the test argument at scripts/manuscript/program_observability_map/analysis_lib.R:538. "
             "The H85 earlier-fraction call (Rho-GTPase / adhesion, q 0.0485) also has tau2 = 0 and Q 0.186 and is 0 "
             "under the truncated variant.")

# ---------------------------------------------------------------- R13-03
e = out['R13-03']
e['changed_in_repair'] = True
e['what_changes'] = (
    "Fig 4E becomes the 774-unit three-cohort arm (364 units): no program meets BH q<0.05, with or without "
    "composition adjustment under either fraction set. 'Sixteen fibrosis and one NAS association ... three fibrosis "
    "and zero NAS' is removed as the main result in the legend, fig4.md and PAPER. The four-cohort arm is a "
    "sensitivity: 16/1, then 3/0 (rebuilt) or 1/0 (earlier, H85) after composition, all needing GSE174478 and the "
    "untruncated Knapp-Hartung variance. The legend now describes the panel as drawn: x is the effect, y is -log10 "
    "BH q, dotted line at q=0.05, no colour scale. Conflict 2: every place that names the 4E arm now names the "
    "364-unit arm (STATUS 6, 585, 608, 814, 888; ROADMAP 299; docs/README 6; fig4 README 8 and 37; INDEX 38; "
    "captions 8), the vocabulary null is quoted for the 364-unit arm (1 vs 1.24, P=0.473; NAS 0 vs 0.19, P=1) with "
    "the 457-record values as sensitivity, and intro 97 gives the ECM/IGFBP7 bulk fibrosis association as q=0.214 "
    "(0.421 adjusted) in the 364 units. RESULTS 1727 and STATUS 608 record that review round 3 replaced the "
    "2026-10-13 choice.")
e['old_value'] = ('16 fibrosis / 1 NAS at BH q<0.05; 3 / 0 after composition (457 records incl. 93 GSE174478); '
                  'vocabulary null 20 vs 26.9 (P=0.871) quoted as the Figure 4E arm; ECM/IGFBP7 bulk q=0.045, '
                  '0.230 adjusted')
e['new_value'] = ('0 / 0 at BH q<0.05 (smallest q 0.079 / 0.159); 0 / 0 after composition with either fraction set '
                  '(364 units, 774 substrate); vocabulary null 1 vs 1.24 (P=0.473), NAS 0 vs 0.19 (P=1); '
                  'ECM/IGFBP7 bulk q=0.214, 0.421 adjusted')
e['source_tables'] = e['source_tables'] + [
    ROOT + 'figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/stromal_ecm_ductular_by_arm.tsv',
    ROOT + 'figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/source_tables/R10-01_vocabulary_null_summary_all_arms.tsv',
    ROOT + 'Analysis/Hypotheses/H85-20261015T124643Z/RESULT_H85.md']
e['commands'] = e['commands'] + (' Repair: read stromal_ecm_ductular_by_arm.tsv, the F4abcd vocabulary-null summary '
                                 '(ws/F4abcd/jobs/7060, 200 draws, seed 20260812, p=(1+#null>=obs)/201) and RESULT_H85.md; '
                                 'pdftotext and PNG of the first-round 4E candidate for the axis encoding.')
e['note'] = (
    'Install as one set only after the user approves reversing the 2026-10-13 four-cohort choice (RESULTS:1727, '
    'STATUS:608). The edits are written in the post-decision form so that the legend, drafts and status documents '
    'agree once installed. The fig4.md:193, RESULTS:1635 and RESULTS:2146 edits supersede F4abcd R10-01 on the same '
    'lines, and the intro 97 edit supersedes F4abcd R04-10 (conflict 2); F4abcd R04-01 edits the next sentence of '
    'intro 99 and does not overlap. The vocabulary-null counts use unweighted scores and a 228-test family, which '
    'differs from the weighted 224-test Figure 4E family; the text says so.')

# ---------------------------------------------------------------- R13-05
e = out['R13-05']
e['changed_in_repair'] = True
e['what_changes'] = (
    "No panel number changes. The abstract, fig3.md, fig4.md and the Discussion stop saying the continuum 'resolves "
    "coordinated' within-stage variation. The abstract says scores 'still varied with the continuum' (no ambiguous "
    "'it'); fig3.md names the roster (367 disease units, 94 of 116 estimable programs unadjusted, 53 to 76 after "
    "adjustment); fig4.md:126 names the 367 units and their meta coefficients (0.760, 0.707) before the "
    "composition-adjusted ranges; fig4.md:149 says the 43 systems behave like about 14 independent systems with "
    "PC1 at 41% (466-unit roster) instead of inferring one shared axis from 14; the Discussion uses 'tracks the "
    "continuum'. Discussion 207 follows joint-rule row 5 (H52 untestable, H54 supports): 'shared yet "
    "cohort-variable' replaces 'coordinated', and the H54 replication in GSE281797 and GSE225740 is stated with "
    "its two limits.")
e['new_value'] = ('367 disease units: 94/116 programs unadjusted, 53-76 adjusted (A rebuilt 67, A earlier 65); '
                  'Stromal ECM 0.760 to 0.203-0.479; Ductular 0.707 to 0.60-0.73; GSE162694 Stromal ECM 0.172, q 0.233; '
                  '0/43 connected null; 14 effective systems, PC1 41.07% (466 units); H54 replicated in A91 (GSE281797) '
                  'and B87 (GSE225740)')
e['source_tables'] = e['source_tables'] + [
    ROOT + 'figures/candidates/continuum-rare-type-arms-20261012T020132Z/arms/meta_effects_all_arms.tsv',
    ROOT + 'figures/candidates/continuum-composition-774-20261012T012137Z/basis_774/continuum/continuum_meta_effects.tsv',
    ROOT + 'Analysis/Hypotheses/H54-20261015T111917Z/RESULT_H54.md']
e['commands'] = e['commands'] + (' Repair: program call counts by arm from meta_effects_all_arms.tsv (U 94/116; '
                                 'range 53-76 over 11 fits); continuum_meta_effects.tsv disease_only rows for the two '
                                 'focal programs (fixed_projection 0.760/0.707 HC3, signature_pc1 0.762/0.699); '
                                 'system_redundancy.tsv; RESULT_H54.md.')
e['note'] = ('Discussion 207-208 is now edited (claim_level) under joint-rule row 5, with both H54 cohorts named; '
             'the title clause is removed elsewhere and is not touched here. The fig4.md:126 and abstract edits also '
             'carry the R13-11 wording (genetic-cluster-level, samples). Cohort-PC1 concordance on 367 units: '
             'signature_pc1 0.762 (Stromal ECM) and 0.699 (Ductular injury).')

# ---------------------------------------------------------------- R13-11
e = out['R13-11']
e['changed_in_repair'] = True
e['what_changes'] = e['what_changes'] + (" Repair: the fig4.md:1 section header 'resolve cell-state and regulatory "
                                         "context' becomes 'Donor-level single-cell programs and their bulk-tissue "
                                         "context', matching the fig1.md sentence.")
e['new_value'] = e['new_value'] + "; fig4.md header 'Donor-level single-cell programs and their bulk-tissue context'"
e['note'] = (e['note'] + ' Conflict 24 (fig1.md:44, F4abcd R13-04 on the same old text) is not owned by F4efgh; its '
             'proposed merged sentence keeps the bulk-tissue context and supplementary TF wording used here.')

entries = [out[f] for f in ['R05-02', 'R05-04', 'R05-05', 'R05-06', 'R05-07', 'R05-08', 'R05-09',
                            'R12-01', 'R13-03', 'R13-05', 'R13-11']]
for e in entries:
    for k in ('what_changes', 'note', 'new_value', 'reproduced'):
        if k in e:
            assert '—' not in e[k], (e['fid'], k)

cfr = [
    dict(fid='R05-02', value='Fig 4G and 4H on different rosters (S01 0.712 in 4G vs 0.665 in 4H)',
         who_was_right='checker', final_value='4G re-rendered on the 367 disease units; S01 0.665 in both panels',
         evidence=REP + 'fig4g_node_checks.tsv: system betas equal the 4H source to 8.3e-17 (disease_only) and 9.0e-17 '
                  '(as_published); as_published nodes reproduce installed 4G to 1.3e-15; PDF ' + G4_PDF +
                  ' (sha256 0e936337b580...); legend (G)+(H) and fig4.md:156 name the 367-unit roster.'),
    dict(fid='R05-02', value="PAPER:32 'Figure 4H plots the 367 disease units'",
         who_was_right='checker', final_value="'Figure 4H plots 43 system coefficients estimated in the 367 disease units'",
         evidence='Panel source fig4h_disease_only_source.tsv has 43 rows (systems), not units.'),
    dict(fid='R05-04', value='GSE162694 earlier-fraction Stromal ECM upper bound 0.639',
         who_was_right='checker', final_value='0.638 (0.63850)',
         evidence='meta_effects_all_arms.tsv and checker job 9073: 0.4334 [0.2283, 0.6385].'),
    dict(fid='R05-06', value="'most S4U asterisks mark systems whose call does not survive adjustment'",
         who_was_right='checker',
         final_value='4/18 called in all 11 fits, 14/18 lose the call in at least one, 9/18 called in 5 or fewer',
         evidence='Per-system call counts over 11 estimable fits: 0:1, 2:2, 3:2, 4:1, 5:3, 6:1, 10:4, 11:4 '
                  '(meta_effects_all_arms.tsv joined with fig4h_source.tsv supported==True).'),
    dict(fid='R05-06', value="'With about 14 effectively independent systems, these counts describe one shared axis'",
         who_was_right='checker',
         final_value='Counts reported as part of a diffuse shared axis because system scores are strongly correlated '
                     'and 0/43 exceed the connected-region null; 14 effective systems and PC1 41% stated as facts',
         evidence='system_redundancy.tsv (466 units): 14 effective systems pooled, PC1 41.07%; s4w gate connected 0/43.'),
    dict(fid='R12-01', value="'both 774-unit arms give 0 calls'; z-test 'not the prespecified test'",
         who_was_right='checker',
         final_value="'both arms give 0 calls'; 'the Figure 4E code uses Knapp-Hartung'",
         evidence='The four-cohort arm adds 93 GSE174478 records outside the 774 units; test="knha" at '
                  'scripts/manuscript/program_observability_map/analysis_lib.R:538; no prespecification text found.'),
    dict(fid='R13-03', value="Legend (E) 'color shows BH q over 224 tests'",
         who_was_right='checker',
         final_value='x axis effect, y axis -log10 BH q over 224 tests, dotted line at q = 0.05',
         evidence='pdftotext and PNG of ' + CAND + 'fig4e/panels/fig4e_bulk_tissue_state_transport.pdf: one gray, '
                  '112 points per facet, dotted q line, no colour key.'),
    dict(fid='R13-03', value="RESULTS:1727 and STATUS:609 say 'user decision pending' while legend, fig4.md and PAPER present the 364-unit arm",
         who_was_right='checker',
         final_value='All edits written in the post-decision form, to be installed as one set after user approval',
         evidence='Conflict 2 resolution; RESULTS:1727 and STATUS:608 edits now record that round 3 replaced the 2026-10-13 choice.'),
    dict(fid='R13-05', value="Abstract 'varied with it along one shared axis'",
         who_was_right='checker', final_value="'still varied with the continuum'",
         evidence='Wording only; antecedent of it was ambiguous between stage and continuum.'),
    dict(fid='R13-05', value="fig4.md:126 unlabeled 367-unit values beside 435-unit models",
         who_was_right='checker',
         final_value='Sentence opens with the 367 disease units and gives their meta coefficients 0.760 and 0.707',
         evidence='continuum_meta_effects.tsv disease_only, n=367: fixed_projection HC3 0.760 (Stromal ECM), 0.707 '
                  '(Ductular); signature_pc1 0.762, 0.699.'),
    dict(fid='R13-05', value="fig3.md '94 programs called without adjustment' and 'the 117 programs'",
         who_was_right='checker',
         final_value='367 disease units, 94 of 116 estimable programs unadjusted, 53 to 76 after adjustment',
         evidence='meta_effects_all_arms.tsv: U arm 94/116 estimable; adjusted fits 53 to 76 (A rebuilt 67, A earlier 65).'),
    dict(fid='R13-11', value="fig4.md:1 header 'resolve cell-state and regulatory context' left unedited",
         who_was_right='checker', final_value="'## Donor-level single-cell programs and their bulk-tissue context'",
         evidence='grep of docs/manuscript/draft/fig4.md line 1; added as an R13-11 edit.'),
]

conf = [
    dict(index=2, resolution=(
        'Figure 4E is the 364-unit three-cohort arm everywhere it is named, written in post-decision form under '
        'R13-03: legend 248, fig4.md 5/113/193, PAPER 444, intro 97, RESULTS 1635/1713/1727/2146, STATUS 6/585/608/814/888, '
        'docs/README 6, ROADMAP 299, fig4 README 8 and 37, INDEX 38, captions 8. fig4.md 193 reads 1 vs 1.24 (P=0.473), '
        'NAS 0 vs 0.19 (P=1), with the 457-record arm as sensitivity (20 vs 26.9, P=0.871; 2 vs 2.57, P=0.687); intro 97 '
        'gives q=0.214 (0.421 adjusted) in the 364 units and q=0.045 only with 93 GSE174478 records. These supersede '
        'F4abcd R10-01 (fig4.md 193, RESULTS 1635, 2146) and R04-10 (intro 97). The H85 earlier-fraction count (1/0) is '
        'included. This reverses the user choice of 2026-10-13, so the set is installed only after user approval.')),
    dict(index=13, resolution=(
        'fig4 README 8-12 now reads rho=0.63 and 32/43 on the 367 disease units with no outline (row 40 likewise); '
        'ROADMAP 282 gives 0.63 and 32/43 on 367 and says the 435-unit maxT and competitive counts include the 68 '
        'GSE213621 controls; STATUS 303, INDEX 38 and captions 8 updated to the same values. Edits listed under R05-02 '
        '(ROADMAP 282, README 39/40, STATUS 303) and R13-03 (README 8, INDEX, captions).')),
    dict(index=14, resolution=(
        'Legend 203 adds that the 435 units include the 68 GSE213621 controls with no recorded histology placed at '
        'stage 0 and that 4G/4H use 367; legend 317 says the matched S3 comparisons use those 435 units; the RESULTS '
        'layer table gains a 367-unit effect rho column (Hotspot 0.629, Hallmark 0.662, KEGG 0.675, Reactome 0.612, '
        'GO BP 0.671, GO CC 0.547, GO MF 0.523; counts not recomputed). Figure 4G was re-rendered on 367 units so 4G '
        'and 4H share one roster: ' + G4_PDF + '. All listed under R05-02.')),
]

res = dict(ws='F4efgh', entries=entries, check_failures_resolved=cfr, conflicts_resolved=conf)
json.dump(res, open(HERE + '/repair_output_F4efgh.json', 'w'), indent=1)
print('entries', len(entries), 'changed', sum(e['changed_in_repair'] for e in entries),
      'text_edits', sum(len(e.get('text_edits', [])) for e in entries))
