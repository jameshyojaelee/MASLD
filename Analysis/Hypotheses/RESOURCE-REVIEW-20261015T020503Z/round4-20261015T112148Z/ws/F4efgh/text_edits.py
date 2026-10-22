"""Drafted text edits for workstream F4efgh. Validates that each old_text occurs exactly once
in its file and records the 1-based start line. Writes text_edits.json. Edits nothing."""
import json
from pathlib import Path

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
LEG = 'docs/manuscript/05_figure_legends.md'
F4 = 'docs/manuscript/draft/fig4.md'
E = []


def ed(fid, f, old, new, claim=False):
    E.append(dict(fid=fid, file=f, old_text=old, new_text=new, claim_level=claim))


# ---------------- R13-03 / R12-01: Figure 4E ----------------
ed('R13-03', LEG, """**(E)** Weighted program-score linear associations with recorded fibrosis
category and NAS in the four-cohort transport arm: 457 records comprising
364 selected units from GSE130970, GSE135251 and GSE162694 plus 93 source
records from GSE174478. This arm is separate from the five-cohort 774-unit
disease contrast. Effects are score SD per stage-code unit or NAS point;
color shows BH q over 224 tests (112 eligible programs × two axes).
Sixteen fibrosis and one NAS association meet BH q <0.05; after composition
adjustment, three fibrosis and zero NAS associations remain. GSE174478's
source-record identity qualification remains separate from the selected
genetic clusters. Linear category coding does not establish nonlinear or
longitudinal remodeling.""", """**(E)** Weighted program-score linear associations with recorded fibrosis
category and NAS on the 774-unit substrate: 364 selected genetic-identity
clusters from GSE130970, GSE135251 and GSE162694, the three Resource cohorts
that record both axes. Effects are score SD per stage-code unit or NAS point
from a three-cohort random-effects meta-analysis (REML, Knapp-Hartung); color
shows BH q over 224 tests (112 eligible programs × two axes). No fibrosis or
NAS association meets BH q <0.05 (smallest q 0.079 and 0.159), with or without
composition adjustment. Adding 93 GSE174478 source records, which lie outside
the 774-unit substrate and are not identity-resolved, gives 16 fibrosis and
one NAS association at q <0.05 (three and zero after composition
adjustment); all 20 have between-cohort Q below its expected value, and none
remains when the Knapp-Hartung variance is not allowed to fall below the REML
random-effects variance. Linear category coding does not establish nonlinear or
longitudinal remodeling.""")
ed('R13-03', F4, """Figure 4E separately
> includes 93 source records from GSE174478 outside that substrate.""", """Figure 4E uses the
> 774-unit substrate only; 93 GSE174478 source records enter a sensitivity arm.""")
ed('R13-03', F4, """Figure 4E shows linear fibrosis-category and NAS associations for all 112
eligible weighted programs, with BH correction over 224 tests. The selected
four-cohort arm uses 364 units from three corrected bulk cohorts plus 93
GSE174478 source records, 457 in total. Sixteen fibrosis and one NAS association
meet BH q<0.05; after composition adjustment, three fibrosis and no NAS
association remain. The arm differs from the five-cohort pooled disease
contrast and retains GSE174478's source-record identity limitation. Linear
category coding does not summarize nonlinear remodeling or within-person
change.""", """Figure 4E shows linear fibrosis-category and NAS associations for all 112
eligible weighted programs, with BH correction over 224 tests, in the 364
selected genetic-identity clusters of the three 774-unit cohorts that record
both fibrosis and NAS. No program meets BH q<0.05 on either axis (smallest q
0.079 for fibrosis and 0.159 for NAS), with or without composition adjustment.
Adding 93 GSE174478 source records, which are outside the 774-unit substrate
and not identity-resolved, raises the meta-analysis from three to four cohorts
and gives 16 fibrosis and one NAS call (three and none after composition
adjustment). All 20 calls depend on the untruncated Knapp-Hartung variance:
each has between-cohort Q below its expected value, and none remains when that
variance is floored at the REML value. Figure 4E therefore shows no program
association with linear fibrosis category or NAS that holds across these
analysis choices. Linear category coding does not summarize nonlinear
remodeling or within-person change.""")
ed('R13-03', 'docs/PAPER.md', """The selected four-cohort arm uses GSE130970, GSE135251 and GSE162694
(364 selected genetic-identity units) plus GSE174478 (93 source records),
457 in total. GSE174478 is outside the 774-unit pooled disease fit and
retains a separate source-record identity limitation. Sixteen fibrosis and
one NAS association meet BH q<0.05; after composition adjustment, three
fibrosis and zero NAS associations remain. The three-cohort sensitivity
has no call.""", """Figure 4E uses GSE130970, GSE135251 and GSE162694 (364 selected
genetic-identity units), the three 774-unit cohorts that record fibrosis and
NAS. No fibrosis or NAS association meets BH q<0.05 (smallest q 0.079 and
0.159), with or without composition adjustment. A sensitivity arm that adds
93 GSE174478 source records (outside the 774-unit fit, not identity-resolved)
gives 16 fibrosis and one NAS call, three and zero after composition
adjustment; all 20 depend on the untruncated Knapp-Hartung variance (Q below
its expectation; none under the truncated variant) and are not reported as
findings.""")
ed('R13-03', 'docs/RESULTS.md', """The 2026-10-13 Figure 4 candidate uses the four-cohort arm (774 units plus
GSE174478, named in the legend), which the main session recommended and the
user told it to proceed with. The three-cohort version is kept beside it.""", """The 2026-10-13 Figure 4 candidate used the four-cohort arm (774 units plus
GSE174478, named in the legend), which the main session recommended and the
user told it to proceed with; both arms had been computed on 2026-10-12 before
that choice. Review round 3 (R13-03, R12-01) recommends the three-cohort
774-unit arm as Figure 4E (0 fibrosis and 0 NAS calls) and reports the
four-cohort arm as a sensitivity whose calls need GSE174478 and untruncated
Knapp-Hartung (candidate
`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4e/`); user
decision pending.""")
ed('R13-03', 'docs/STATUS.md', """recommendations: Figure 4E is the four-cohort arm (774 units plus GSE174478,
named in the legend), and Figure 4H""", """recommendations: Figure 4E is the four-cohort arm (774 units plus GSE174478,
named in the legend; review round 3 recommends the three-cohort 774-unit arm,
0 and 0 calls, because all 20 four-cohort calls need GSE174478 and the
untruncated Knapp-Hartung variance; user decision pending), and Figure 4H""")
ed('R12-01', 'docs/RESULTS.md', """The loss of calls in the three-cohort arm follows from the cohort count: the
old records without GSE174478 (376 records, three cohorts) give 1 fibrosis and
0 NAS calls.""", """The four-cohort calls are not explained by the cohort count alone. All 20
(16 fibrosis, 1 NAS, 3 composition-adjusted fibrosis) have tau2 = 0 and
between-cohort Q of 0.09 to 1.61, below the expected value of 3, so their
untruncated Knapp-Hartung SE is 0.17 to 0.73 of the REML SE. With the
Knapp-Hartung variance floored at the REML value (metafor test="adhoc"), both
774-unit arms give 0 calls in every model (smallest q 0.193 in the four-cohort
arm). A REML z-test gives 47/28/21/10 calls in the four-cohort arm and
37/14/13/11 in the three-cohort arm; with three or four cohorts it ignores the
uncertainty in tau2 and was not the prespecified test
(`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4e/source_tables/fig4e_calls_by_arm_and_hk_variant.tsv`).
The old records without GSE174478 (376 records, three cohorts) give 1 fibrosis
and 0 NAS calls.""")

# ---------------- R05-02: Figure 4H stage axis ----------------
ed('R05-02', LEG, """are unchanged. **(H)** Fibrosis-stage coefficient (score SD per cohort
stage-code SD) versus stage-adjusted continuum coefficient (score SD per
continuum SD) for all 43 systems, using 435 stage-complete units (109
GSE162694, 326 GSE213621). Spearman rho is 0.7265 and 37/43 effects share
direction; both summaries are descriptive. One outline identifies the system
exceeding the uniform size-matched competitive null on both co-primary axes;
none exceeds the connected-region null. The 18 all-four strict maxT calls are
a different test and are not outlined in the main panel.""", """are unchanged. **(H)** Fibrosis-stage coefficient (score SD per cohort
stage-code SD) versus stage-adjusted continuum coefficient (score SD per
continuum SD) for all 43 systems, using the 367 stage-complete disease units
(109 GSE162694, 258 GSE213621). The 68 GSE213621 controls have no recorded
histology and are excluded from both axes. Spearman rho is 0.63 and 32/43
effects share direction; both summaries are descriptive. The system maxT and
competitive-null tests were run on 435 units that include those controls and
are not marked in this panel.""")
ed('R05-02', F4, """continuum effects have Spearman rho=0.7265 and 37 share direction. The
stage-complete evaluation roster contains 435 units, 109 GSE162694 and 326
GSE213621.""", """continuum effects have Spearman rho=0.63 and 32 share direction. The
stage-complete evaluation roster contains 367 disease units, 109 GSE162694 and
258 GSE213621; the 68 GSE213621 controls have no recorded histology and are
excluded.""")
ed('R05-02', 'docs/manuscript/METHODS.md', """No step producing a
Figure 4G or 4H value reads a control label or disease-versus-control result;
fixed system memberships and outcome-blind coordinates remain unchanged.""", """In the five-cohort
manifest the GSE213621 stage code places its 68 controls, which have no
recorded histology, at stage 0. Figure 4H stage and continuum coefficients are
therefore estimated in disease units only (109 GSE162694 and 258 GSE213621).
The Figure 4G overlay and the system maxT and competitive-null tests were run
on 435 units in which those controls form their own stage category. Fixed
system memberships and outcome-blind coordinates remain unchanged.""")
ed('R05-02', 'docs/RESULTS.md', """rho=0.658 and 36/43 share direction (844 records; 0.7265 and 37/43 on the 774
units, see Main Figure 4H below)""", """rho=0.658 and 36/43 share direction (844 records; 0.7265 and 37/43 on the 774
units with the 68 GSE213621 controls coded as stage 0; 0.631 and 32/43 on the
367 disease units; see Main Figure 4H below)""")
ed('R05-02', 'docs/RESULTS.md', """The stage-adjusted models use 435 units
(109 GSE162694 and 326 GSE213621;""", """The stage-adjusted models use 435 units
(109 GSE162694 and 326 GSE213621, the latter including 68 controls that the
stage code places at stage 0. With those controls removed from both axes the
stage-versus-continuum effect rho values below become 0.629 for Hotspot
programs (81/116 same direction), 0.662 Hallmark, 0.675 KEGG, 0.612 Reactome,
0.671 GO BP, 0.547 GO CC and 0.523 GO MF; the stage-supported counts have not
been recomputed;""")
ed('R05-02', 'docs/RESULTS.md', """Stage-adjusted n=435 (109 +326); these are""", """Stage-adjusted n=435 (109 +326, including the 68 GSE213621 controls, which
carry no recorded histology); these are""")
ed('R05-02', 'docs/RESULTS.md', """and no step that produces a Figure 4G or 4H value reads a control label or a disease-versus-control result; these values stand on the corrected inputs""", """but the GSE213621 stage code places its 68 controls at stage 0, so the installed Figure 4H stage axis reads the control label and the continuum and system models include those controls as their own stage category (review round 3, R05-02); the other values stand on the corrected inputs""")
ed('R05-02', 'docs/RESULTS.md', """On the 774 units the panel source gives Spearman rho 0.7265
(printed as 0.73) and 37/43 in the same direction""", """On the 774 units the panel source gives Spearman rho 0.7265
(printed as 0.73) and 37/43 in the same direction; with the 68 GSE213621
controls removed from both axes (367 disease units) the values are 0.631 and
32/43 (candidate `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`)""")
ed('R05-02', 'docs/PAPER.md', """Figure 4H carries one outline, on the single system exceeding
the uniform size-matched null (0 connected, 1 uniform); the 18 strict maxT
systems are not outlined.""", """Figure 4H plots the 367 disease units and carries no outline;
on the 435-unit roster that includes the GSE213621 controls, one system
exceeds the uniform size-matched null (0 connected, 1 uniform), and the 18
strict maxT systems are not outlined.""", True)
ed('R05-02', 'docs/PAPER.md', """Its single outline marks the one system exceeding 5,000 size-matched
draws from the same membership geometry on both co-primary axes. On the 774
units, Spearman rho is 0.7265 and 37/43 effects share direction; both remain
descriptive.""", """Both coefficients are estimated in the 367 stage-complete disease units;
the 68 GSE213621 controls have no recorded histology. Spearman rho is 0.63
and 32/43 effects share direction (0.7265 and 37/43 when the controls entered
as stage 0); both remain descriptive. The size-matched competitive null was
run on the 435-unit roster with the controls, so the panel carries no outline
until it is rerun on disease units.""")
ed('R05-02', 'docs/PAPER.md', """Those models
   are complete-case in recorded stage and therefore use 109 GSE162694 and 326
   GSE213621 units, 435 in total,""", """Those models
   use 109 GSE162694 and 326 GSE213621 units, 435 in total, where the 326
   include 68 controls with no recorded histology that the stage code places at
   stage 0,""")
ed('R05-02', 'docs/STATUS.md', """stage and continuum system effects have Spearman rho
0.73 with 37 of 43 the same sign""", """stage and continuum system effects have Spearman rho
0.73 with 37 of 43 the same sign (0.63 and 32 of 43 on the 367 disease units
once the 68 GSE213621 controls coded as stage 0 are removed; review round 3,
R05-02)""")
ed('R05-02', 'docs/STATUS.md', """Main 4H has one uniform-null outline, rho 0.7265 and 37/43 same direction.""", """Main 4H has one uniform-null outline, rho 0.7265 and 37/43 same direction (the 326 GSE213621 units include 68 controls coded as stage 0; the disease-only candidate has no outline, rho 0.631 and 32/43).""")

# ---------------- R05-04 / R05-05: S3-S4 focal continuum ----------------
ed('R05-05', LEG, """The focal continuum forest shows disease-only stage/sex-adjusted and
composition-adjusted effects on the 774-unit disease-only evaluation roster
(109 GSE162694 and 258 GSE213621 units). The rebuilt 16-type ILR adjustment
reduces the Stromal ECM meta-coefficient from 0.760 to 0.294 (95% CI
0.154–0.434) and Ductular injury from 0.707 to 0.601 (0.500–0.701).
Stromal ECM is not individually supported after adjustment in GSE162694
(0.172, −0.032–0.377; 117-program-family BH q=0.233). Adjustment depends
on the deconvolution input; no composition-adjusted system is highlighted.""", """The focal continuum forest shows stage/sex-adjusted and composition-adjusted
effects on the 774-unit disease-only evaluation roster (109 GSE162694 and 258
GSE213621 units), with the composition adjustment shown for both cell-type
fraction sets. Ductular injury holds after composition adjustment (meta
coefficient 0.60 to 0.73 across codings and fraction sets; unadjusted 0.707).
Stromal ECM is reduced by composition adjustment by an amount that depends on
the deconvolution input (0.203 to 0.479; unadjusted 0.760). In GSE162694
alone, adjusted Stromal ECM is 0.097 to 0.217 with the rebuilt fractions and
0.344 to 0.433 with the earlier fractions. No composition-adjusted program or
system is named as a finding.""")
ed('R05-05', 'docs/STATUS.md', """composition-adjusted 0.294 [0.154, 0.434], BH q 1.4e-4 (0.726 and 0.488 on the
844 records); Ductular injury 0.707 and 0.601 (0.677 and 0.639). Adjusted
Stromal ECM in GSE162694 alone is 0.172 [-0.032, 0.377], BH q 0.233, no longer
supported there.""", """composition-adjusted 0.294 [0.154, 0.434], BH q 1.4e-4 in arm A with the rebuilt
fractions (0.203 to 0.479 across arms and fraction sets; 0.726 and 0.488 on the
844 records); Ductular injury 0.707 and 0.601 (0.60 to 0.73 across arms; 0.677
and 0.639 on the 844 records). Adjusted Stromal ECM in GSE162694 alone is 0.172
[-0.032, 0.377], BH q 0.233, in arm A with the rebuilt fractions (0.097 to 0.217
across arms); with the earlier fractions it is 0.344 to 0.433.""")

# ---------------- R05-06 / R05-07: strict maxT 18 ----------------
ed('R05-06', LEG, """Figure 4H uses one uniform-null outline, as described above. The supplementary
inference registry retains 33 two-axis BH/direction calls and 18 all-four
strict maxT calls; the jointly powered cross-cohort maxT test gives 33/43.""", """The supplementary inference registry, fitted on 435 units in which the 68
GSE213621 controls form their own stage category, retains 33 two-axis
BH/direction calls and 18 strict maxT calls (maxT P<0.05 in both cohorts on
both co-primary axes); the jointly powered cross-cohort maxT test gives 33/43.
The two co-primary axes give system coefficients correlated at r=0.999, so the
strict criterion is in effect a two-cohort conjunction. The 18 strict calls
come from models without composition terms: 4 of the 18 are called (BH q<0.05
within 43 systems) in all 11 estimable composition-adjusted fits and 1 in
none, so most S4U asterisks mark systems whose call does not survive
adjustment for cell-type composition.""")
ed('R05-06', 'docs/PAPER.md', """S4U retains
   all 43 systems and marks only those 18 strict calls with a neutral asterisk;""", """S4U retains
   all 43 systems and marks only those 18 strict calls with a neutral asterisk
   (4 of the 18 are called in all 11 estimable composition-adjusted fits, and
   1 in none);""")
ed('R05-06', F4, """Figure 4H compares the stage coefficient in score SD per cohort SD of the
recorded stage code with the stage/sex-adjusted continuum coefficient in score
SD per continuum SD. These estimate different quantities. Its single outline
marks the system exceeding the uniform size-matched competitive null on both
co-primary axes; no system exceeds the connected-region null on both axes.
The separate 18 all-four strict maxT calls remain in the supplement and do not
establish competitive specificity. The jointly powered cross-cohort maxT
test gives 33 of 43 systems on the corrected counts.""", """Figure 4H compares the fibrosis-stage coefficient in score SD per cohort SD
of the recorded stage code with the stage/sex-adjusted continuum coefficient in
score SD per continuum SD, both estimated in the 367 disease units. These
estimate different quantities. The system maxT and competitive-null tests were
run on 435 units that included the 68 GSE213621 controls. On that roster no
system exceeds the connected-region null on both axes, one (S01) exceeds the
permissive uniform size-matched null, and 18 pass maxT in both cohorts on both
co-primary axes. The two axes give system coefficients correlated at r=0.999,
so this is in effect a two-cohort criterion, and only 4 of the 18 are called
in all 11 estimable composition-adjusted fits. The jointly powered
cross-cohort maxT test gives 33 of 43. With about 14 effectively independent
systems, these counts describe one shared axis, not system-specific
associations.""")
ed('R05-07', 'docs/RESULTS.md', """The strict criterion is the conjunction of four
independently maxT-corrected axis-by-cohort tests, and GSE162694 contributes 109""", """The strict criterion is the conjunction of four
maxT-corrected axis-by-cohort tests, but the two co-primary axes give system
coefficients correlated at r=0.9993 (0.9997 and 0.9989 within cohort), and
requiring both axes removes none of the 18 (one axis in both cohorts gives 18
on either axis), so it is in effect two cohort tests. GSE162694 contributes 109""")

# ---------------- R05-09 ----------------
ed('R05-09', 'docs/RESULTS.md', """donor-first and plotted pooled-window estimands agree at Spearman rho=0.997.""", """donor-first and plotted pooled-window estimands agreed at Spearman rho=0.997
on the 844 records. No like-for-like 774-unit reconciliation was run: the
774-unit donor-first windows against the frozen 844-record feature-first
windows give Spearman 0.849 (Pearson 0.922) for the equal-weight five-cohort
pool.""")

# ---------------- R13-05 (and R13-11 wording folded in) ----------------
ed('R13-05', 'docs/manuscript/draft/abstract.md', """An externally anchored molecular continuum
complemented recorded fibrosis stage, preserving the major metabolic,
inflammatory, ductular, and extracellular-matrix patterns while resolving
coordinated transcriptional variation among participants assigned the same
stage.""", """An externally anchored molecular continuum
complemented recorded fibrosis stage, preserving the major metabolic,
inflammatory, ductular, and extracellular-matrix patterns. Among samples
sharing a recorded stage, many program and pathway scores varied with it along
one shared axis; which programs carry that variation depends on adjustment for
cell-type composition.""", True)
ed('R13-05', 'docs/manuscript/draft/fig3.md', """The
continuum recovered the broad stage-associated direction of transcript,
pathway, NMF, and Hotspot remodeling while retaining substantial stage-adjusted
signal, particularly for multicellular programs.""", """The
continuum recovered the broad stage-associated direction of transcript,
pathway, NMF, and Hotspot remodeling while retaining stage-adjusted signal;
for the 117 programs that signal depends on cell-type composition adjustment
and on the deconvolution input (94 programs called without adjustment, 67
after adjustment with the rebuilt fractions).""")
ed('R13-05', F4, """In GSE162694 and GSE213621, both Stromal ECM and Ductular injury increase
along the fixed external projection in stage- and sex-adjusted donor-level
models, with concordant effects for the cohort-PC1 score.""", """In GSE162694 and GSE213621, both Stromal ECM and Ductular injury increase
along the fixed external projection in stage- and sex-adjusted
genetic-cluster-level models, with concordant effects for the cohort-PC1
score. After cell-type composition adjustment, Ductular injury holds (0.60 to
0.73 across codings and fraction sets; unadjusted 0.707), whereas Stromal ECM
falls to 0.203 to 0.479 (unadjusted 0.760) depending on the deconvolution
input; in GSE162694 alone it is not supported with the rebuilt fractions
(0.172 [-0.032, 0.377], BH q 0.233).""")
ed('R13-05', F4, """histology
anchors severity, whereas the continuum resolves coordinated molecular
heterogeneity among participants sharing a recorded histologic category.""", """histology
anchors severity, whereas the continuum describes a shared axis of molecular
variation among samples sharing a recorded histologic category.""")
ed('R13-05', F4, """auditable groupings so that remodeling within a recorded stage can be described
as coordinated rather than gene by gene.""", """auditable groupings so that remodeling within a recorded stage can be described
by feature group rather than gene by gene. No system exceeds the
connected-region null on both axes (0/43), and the 43 systems have about 14
effectively independent components, so system counts describe one diffuse
shared axis rather than 43 separate findings.""")
ed('R13-05', 'docs/manuscript/draft/fig5_discussion.md', """Its strong alignment with fibrosis-based pathway, NMF, and Hotspot
results shows that the major biological programs are preserved, while
stage-adjusted associations show that histologic bins compress coordinated
molecular variation.""", """Its strong alignment with fibrosis-based pathway, NMF, and Hotspot
results shows that the major biological programs are preserved, while
stage-adjusted associations show molecular variation within histologic bins
along one shared axis; which programs carry it depends on adjustment for
cell-type composition.""")

# ---------------- R13-11 ----------------
ed('R13-11', 'docs/manuscript/draft/fig1.md', """Figure 4
resolves donor-level multicellular programs and regulatory activity.""", """Figure 4
resolves donor-level single-cell programs and their bulk-tissue context;
transcription-factor activity is supplementary.""")

# ---------------- validate ----------------
bad = []
for e in E:
    text = (ROOT / e['file']).read_text()
    n = text.count(e['old_text'])
    if n != 1:
        bad.append((e['fid'], e['file'], n, e['old_text'][:60]))
        continue
    e['line'] = str(text[:text.index(e['old_text'])].count('\n') + 1)
    assert '—' not in e['new_text'], 'em dash'
print('edits', len(E), 'bad', bad)
assert not bad
out = Path(__file__).with_name('text_edits.json')
out.write_text(json.dumps(E, indent=1, ensure_ascii=False))
