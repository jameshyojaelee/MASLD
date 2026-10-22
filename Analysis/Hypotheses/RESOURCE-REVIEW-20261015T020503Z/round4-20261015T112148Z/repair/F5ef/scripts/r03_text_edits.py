#!/usr/bin/env python3
"""F5ef repair: assemble the repaired text edits and check them against the CURRENT files.

Checks: each old_text occurs exactly once in its file; no new_text contains an em dash or an
evidence tag; no new_text keeps the disputed wording. Writes out/text_edits_repair.json.
Edits not listed as changed are copied unchanged from the first round (repair_inputs/F5ef.json).
"""
import json
import re
from pathlib import Path

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
R4 = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
OUT = R4 / "repair/F5ef/out"
first = {e["fid"]: e for e in json.load(open(R4 / "repair_inputs/F5ef.json"))["entries"]}


def old(fid, i):
    return first[fid]["text_edits"][i]


LEG = "docs/manuscript/05_figure_legends.md"
FIG5 = "docs/manuscript/draft/fig5.md"
METH = "docs/manuscript/METHODS.md"
PAPER = "docs/PAPER.md"
RES = "docs/RESULTS.md"
ABS = "docs/manuscript/draft/abstract.md"

new = {}

new[("R01-03", 0)] = dict(file=LEG, line="408-409", claim_level=False, old_text=old("R01-03", 0)["old_text"], new_text="""\
**(E)** Genetic-chromatin and protein-histology context for the complete
registry of 117 frozen programs; rows are not selected by disease association.
Rows show the 45 programs with at least one member gene in an adopted COLOC r2
primary pair (the 767 pairs of panel D) that places at least half of its
colocalization posterior in open chromatin of the program's own lineage in at
least one snATAC source (28 programs), or with a robust protein-histology
association (26 programs); 9 programs meet both. All 117 rows are in the source
table; display filtering triggers no re-testing. Bars count genes once each;
dark shading marks a gene that is lineage-open through a direct MASLD/liver-fat
pair (RORA, in 3 programs), and light shading a gene that is lineage-open only
through liver-enzyme pairs. Intervals are rank-residual partial Spearman
correlations between the PXD051911 program protein score and a histologic-burden
axis (first principal component of steatosis, ballooning and lobular
inflammation grades and fibrosis stage) in 46 MASLD participants, adjusted for
acquisition batch, age, BMI and sex, with participant-bootstrap 95% CIs.
Freedman-Lane permutation P (9,999) is BH-corrected over the 94 programs with at
least eight measured member proteins and 20% retained weight; "not scored" marks
the two displayed programs below these thresholds (Elastic-fiber ECM, 6
measured proteins; Memory/MAIT, 3). Robust means q<0.05 with the equal-weight
score, the score without its top gene and both acquisition batches agreeing in
sign. Stromal ECM and Ductular injury are absent: neither contains a colocalized
gene and neither can be scored in protein. Spatial organization is in panel F.""")

new[("R01-03", 1)] = dict(file=FIG5, line="39", claim_level=False, old_text=old("R01-03", 1)["old_text"], new_text="""\
gene, and absent overlap is unresolved context rather than a negative test.

Figure 5E places all 117 frozen programs in two assay-native contexts without
selecting them by disease association. Twenty-eight programs contain at least
one gene from an adopted colocalized pair whose posterior lies mainly in open
chromatin of the program's own lineage. In three of them that gene is RORA,
open through a direct MASLD pair; in the other 25 the open genes come through
liver-enzyme pairs only. In liver protein from 46 MASLD participants, 26 of 94
scorable programs are robustly associated with a histologic-burden axis across
participants (BH over 94): 16 negatively, including oxidative-phosphorylation,
xenobiotic and amino-acid degradation (TAT) programs, and 10 positively,
including fibroblast ECM-proteoglycan and TGFB/collagen programs. The same
protein cohort did not contribute to program discovery. Neither spatial program
in Figure 5F contains a colocalized gene or can be scored in protein.""")

new[("R01-03", 2)] = dict(file=METH, line="1093-1094", claim_level=False, old_text=old("R01-03", 2)["old_text"], new_text="""\
concordance in this same cohort, so it is descriptive and is not an independent
validation of the bulk or single-cell programs.

For Figure 5E, single-gene protein groups with at most 50% missing values are
log2-transformed, quantile-normalized, collapsed to one value per symbol by the
median and z-scored per protein across the 58 samples. Each of the 117 frozen
programs is scored per participant as the L1-weighted mean of its measured
member proteins; a program is scorable with at least eight measured members and
at least 20% of its L1 weight retained (94 of 117). Histologic burden is the
first principal component of steatosis, ballooning and lobular inflammation
grades and Kleiner fibrosis stage in the 46 MASLD participants with complete
histology; NAS is excluded because it sums three of these. The association is a
rank-residual partial Spearman correlation adjusted for acquisition batch, age,
BMI and sex, with a batch-stratified Freedman-Lane permutation P (9,999
permutations), a participant bootstrap 95% CI resampled within batch (4,999
draws) and Benjamini-Hochberg correction over the 94 scorable programs. A
program is robust when q<0.05 and the equal-weight score, the score without its
top-weighted gene and both acquisition batches agree in sign. Seeds are 20260923
(permutation) and 30260923 (bootstrap) plus the program's position, from 1 to
94, among the scorable programs in registry order.""")

new[("R01-03", 4)] = dict(file="docs/STATUS.md", line="1076-1077", claim_level=False, old_text=old("R01-03", 4)["old_text"], new_text="""\
- Review the 117-program Figure 5E candidate and the corrected Figure 5F source
  table (`repair/fig5f_spatial_program_calibration_corrected_r2.tsv`) in
  `figures/candidates/review-r3-fixes-20261015T112148Z/F5ef/`; positive Moran
  direction alone is not the main spatial result.""")

new[("R06-10", 1)] = dict(file="docs/STATUS.md", line="306", claim_level=False, old_text=old("R06-10", 1)["old_text"], new_text=(
    "GSE192741, 4 donors: 0.139, exact p 0.0625, the smallest possible with 4 donors, outside the family. These values "
    "adjust for the cell2location hepatocyte abundance (GSE192741, Vu) or a 15-marker hepatocyte proxy (HMSMA), not the "
    "16 abundances of Figure 5F (Stromal ECM excess 0.040 and 0.050), so V4 does not replicate the 5F statistic; Vu "
    "arrays are not donors."))

new[("R06-02", 0)] = dict(file=FIG5, line="41-50", claim_level=False, old_text=old("R06-02", 0)["old_text"], new_text="""\
Figure 5F examines tissue organization for two programs fixed in advance from
the superseded August stage fit; neither passes the corrected 117-program rule.
After adjustment for all 16 cell2location abundances, ECM/IGFBP7 exceeds a
matched null built from genes of equivalent expression, detection, biotype and
correlation with hepatocyte abundance in both sources (excess 0.040 and 0.050),
whereas ductular-injury/BICC1 does not (0.015 and 0.003) and is indeterminate
rather than negative. Residual Moran's I is positive for ductular-injury/BICC1
in 4 of 4 donors and 9 of 10 arrays, so the raw statistic does not separate the
programs, and no prespecified test compared them. The matched-gene null holds
donors fixed. At the donor level, ECM/IGFBP7 excess is positive in 3 of 4
donors (exact sign-flip P=0.19); with four donors no sign-flip test can reach
P<0.05. Disease association, cell-of-origin assignment and tissue architecture
are separate claims requiring separate calibration. Donor joins, array-level
source dependence and the matched-null construction are in Methods.""")

new[("G3-03", 0)] = dict(file=PAPER, line="333-346", claim_level=True, old_text=old("G3-03", 0)["old_text"], new_text="""\
For spatial omics, the finding concerns two programs fixed in advance from the
superseded August stage fit; neither passes the corrected 117-program rule.
**Reproducibly positive spatial autocorrelation is not sufficient evidence that
a program has program-specific tissue architecture.** After residualizing for
all 16 cell2location abundances, both programs had positive aggregate residual
Moran's I in both Visium sources, and residual Moran's I was positive in 3/4 and
4/4 GSE192741 donors and 10/10 and 9/10 Vu arrays for Stromal ECM and Ductular
injury. Stromal ECM exceeded its matched-gene null in both sources (excess
0.040, q=0.001; 0.050, q=0.0006) and Ductular injury did not (0.015, q=0.079;
0.003, q=0.32); no prespecified test compared the two programs. The null
resamples genes matched on expression, detection, biotype and correlation with
hepatocyte abundance, with donors held fixed. At the donor level, Stromal ECM
excess was positive in 3 of 4 donors (exact one-sided sign-flip P=0.19); four
donors cannot reach P<0.05 (smallest attainable P 0.0625). The Vu result remains
source-dependent because arrays cannot be joined to donors. Calibration against
an assay's intrinsic spatial propensity, not the sign or apparent replication of
Moran's I, is therefore required before spatial organization is called
program-specific.""")

new[("G3-03", 2)] = dict(file=RES, line="3150", claim_level=False, old_text=old("G3-03", 2)["old_text"], new_text="""\
sections had residual Moran's I of 0.080 and 0.328 with hepatocyte abundance as
the only composition covariate (0.003 and 0.077 under the all-cell adjustment of
Figure 5F), but they""")

new[("G3-04", 0)] = dict(file=ABS, line="55-59", claim_level=True, old_text=old("G3-04", 0)["old_text"], new_text="""\
physical tissue contexts. In two Visium sources, an ECM/IGFBP7 program
exceeded a matched-gene spatial null after adjustment for cell composition and
a ductular-injury/BICC1 program did not; no prespecified test compared the two
programs, and with four donors the ECM/IGFBP7 result does not reach P<0.05 at
the donor level. Disease association, tissue architecture and cellular origin
are reported as separate evidence.""")

new[("G3-04", 1)] = dict(file=PAPER, line="278-281", claim_level=True, old_text=old("G3-04", 1)["old_text"], new_text="""\
After all-cell adjustment, Stromal ECM
   exceeds its matched-gene spatial null in both Visium sources and Ductular
   injury does not (indeterminate); no prespecified test compared the two
   programs. The null holds donors fixed, and Stromal ECM excess is positive in
   3 of 4 donors. Disease remodeling, tissue architecture, and lineage
   attribution are reported as separate claims.""")

new[("G3-04", 2)] = dict(file=PAPER, line="394", claim_level=True, old_text=old("G3-04", 2)["old_text"], new_text=(
    "| Stromal ECM exceeds a matched-gene spatial null after all-cell adjustment and Ductular injury does not "
    "(indeterminate); no prespecified test compared the two programs, and disease association, tissue "
    "architecture, and cellular origin are reported as separate claims |"))

new[("R12-03", 0)] = dict(file=LEG, line="409-417", claim_level=False, old_text=old("R12-03", 0)["old_text"], new_text="""\
**(F)** At
left, a registered H&E overlay shows the all-cell-adjusted Stromal ECM score in
one Guilliams et al. steatotic section (JBO019, donor H37) selected outcome-blind
as closest to the dataset median eligible spot count; of the five sections, it
has the largest Stromal ECM excess (Results). At right, points show Moran excess
over matched genes for the two prespecified programs in Guilliams et al. and Vu
et al., after residualizing spot scores for all 16 cell2location abundances, log
library size and log detected genes; labels are BH-adjusted empirical q-values
from 9,999 matched-gene draws, BH over the two programs within each source. The
null resamples genes with donors and arrays held fixed, so the q-values are not
donor-level tests; donor-level values are in Results. Stromal ECM exceeds its
null in both sources (supported in Guilliams et al.; source-dependent in Vu et
al. because its arrays lack a public donor key), whereas Ductular injury is
indeterminate in both.""")

new[("R12-07", 0)] = dict(file=RES, line="2994-3017", claim_level=False, old_text=old("R12-07", 0)["old_text"], new_text="""\
The confirmatory family was recorded before execution and contains only the Stromal
ECM and Ductular injury programs, fixed in advance from the superseded August
stage fit; neither passes the corrected 117-program rule. The values below are
those plotted in Figure 5F: spot scores residualized for all 16 cell2location q05
abundances plus log library size and log detected genes, expressed as excess
over the matched-gene null mean, with 9,999 matched-gene draws (seed 42 for
GSE192741, 1042 for Vu) and BH over the two programs within each source. In
GSE192741, Stromal ECM had donor-first excess 0.04050 (matched-null SD 0.01098;
4 of 9,999 draws at or above it, P=0.0005, q=0.001) and is `supported`.
Ductular injury had excess 0.01532 (SD 0.01060, P=q=0.0793) and is
`indeterminate`, not negative.

In Vu, Stromal ECM had excess 0.05034 (SD 0.01033; 2 of 9,999 draws, P=0.0003,
q=0.0006). Its evidence state is `source_dependent` because the ten arrays cannot
be joined to biological donors. The r8 Catalog program-context row for Vu
Stromal ECM reads `supported` and should read `source_dependent`. Ductular
injury had excess 0.002654 (SD 0.008228, P=q=0.3158) and is `indeterminate`.

With the cell2location hepatocyte abundance as the only composition covariate
(the discovery lineage of both programs), the excesses were larger: Stromal ECM
0.1390 in GSE192741 and 0.1471 in Vu, each with no null draw at or above it
(P ≤ 1e-4, the resolution floor of 9,999 draws; q ≤ 0.0002), and Ductular injury
0.02587 (P=q=0.0755) and 0.01922 (P=q=0.1158). This sensitivity changes no
evidence state. The four r8 Catalog spatial program-context rows carry these
hepatocyte-only estimates, not the plotted all-cell values.

The matched-gene null resamples genes; donors and arrays are held fixed. Under
the plotted all-cell model, residual Moran's I was positive for Stromal ECM in
3/4 GSE192741 donors and 10/10 Vu arrays, and for Ductular injury in 4/4 donors
and 9/10 arrays, so neither positive sign nor cross-unit consistency
distinguishes the programs. Stromal ECM excess over each donor's own
matched-gene null was positive in 3 of 4 donors (0.0156, 0.0657, 0.1023 and
-0.0216 for H35 to H38; exact one-sided sign-flip P=0.19, descriptive; the floor
with four donors is 0.0625). Stromal ECM exceeded its matched-gene null in both
sources and Ductular injury did not; no prespecified test compared the two
programs. Descriptively, the per-unit Stromal ECM minus Ductular injury excess
was positive in 3 of 4 GSE192741 donors and 10 of 10 Vu arrays (arrays are not
donors). **Generic spatial autocorrelation can mimic program transportability,
and an assay-matched null is required to separate the two.** This is not a
negative Ductular injury result, and the Vu contrast remains source-dependent
rather than donor-level replication.""")

new[("G3-08", 0)] = dict(file=METH, line="1098-1103", claim_level=False, old_text=old("G3-08", 0)["old_text"], new_text="""\
Each Visium spot receives a program score: the z-scored expression of its
measured fixed member genes, weighted by their original L1 weights rescaled to
sum to one. The score is residualized for all 16 cell2location q05 abundance
factors, log1p library size and log1p detected genes (Figure 5F); a sensitivity
analysis replaces the 16 abundances with the cell2location hepatocyte abundance
alone, the discovery lineage of both programs. Residual Moran's I on six-neighbour
spot graphs built within connected tissue islands is compared with 9,999 random
sets of genes matched on expression, detection rate, biotype and correlation with
the hepatocyte abundance, scored and residualized identically (seed 42 for
GSE192741, 1042 for Vu). The reported excess is the observed statistic minus the
null mean; P=(k+1)/(9,999+1), where k counts draws at or above the observed value,
so the smallest attainable P is 1e-4. Benjamini-Hochberg correction covers the
two prespecified programs within each source; no prespecified test compares the
two programs. The null resamples genes and holds donors and arrays fixed.
Donor-first estimates are used where a donor key exists (five GSE192741 sections
collapse to four donors); an exact sign-flip test of per-donor excess has a
smallest attainable P of 0.0625 with four donors and is reported descriptively.
Array-level results without a donor key are reported as `source_dependent` and
never as donor-level inference.""")

# edits carried unchanged from round 1
unchanged = [("R01-03", 3), ("R11-03", 0), ("R06-02", 2), ("R06-10", 0), ("G3-07", 0)]
# edits dropped: H90 is untestable, so the 'H90 pass' variants do not apply
dropped = [("R06-02", 1), ("G3-03", 1)]

order = ["R01-03", "G3-01", "R11-03", "R06-02", "G3-03", "G3-04", "R12-03", "R06-10", "R06-11", "R12-07", "G3-07", "G3-08"]
result = {}
problems = []
BANNED_NEW = ["neither\nwas selected by a stage-association rule", "neither was selected by a stage-association rule",
              "separated them at the program level", "therefore separates the programs", "0.008227",
              "program's own lineage abundance", "amino-acid metabolism programs", "track a histologic-burden"]
for fid in order:
    edits = []
    n_first = len(first[fid].get("text_edits", []))
    for i in range(n_first):
        if (fid, i) in dropped:
            continue
        if (fid, i) in new:
            e = dict(new[(fid, i)])
            e["changed_in_repair"] = True
        else:
            assert (fid, i) in unchanged, (fid, i)
            e = dict(first[fid]["text_edits"][i])
            e["line"] = e.get("line", "").split(" [")[0]
            e["changed_in_repair"] = False
        text = (P / e["file"]).read_text()
        cnt = text.count(e["old_text"])
        if cnt != 1:
            problems.append(f"{fid}[{i}] old_text count {cnt} in {e['file']}")
        if "—" in e["new_text"]:
            problems.append(f"{fid}[{i}] em dash in new_text")
        if re.search(r"\[(Inference|Unverified|Estimate|Speculation)", e["new_text"]):
            problems.append(f"{fid}[{i}] evidence tag in new_text")
        for b in BANNED_NEW:
            if b in e["new_text"]:
                problems.append(f"{fid}[{i}] keeps disputed wording: {b!r}")
        e["start_line_found"] = text[: text.index(e["old_text"])].count("\n") + 1 if cnt == 1 else None
        edits.append(e)
    result[fid] = edits

json.dump(result, open(OUT / "text_edits_repair.json", "w"), indent=1, ensure_ascii=False)
for fid in order:
    for e in result[fid]:
        print(fid, e["file"], e["line"], "found at", e["start_line_found"], "changed" if e["changed_in_repair"] else "unchanged")
print("dropped:", dropped)
print("PROBLEMS:", problems if problems else "none")
