"""Build RUN/staged/main_misc.json (round-2 writer misc). Checks each old string
occurs exactly once in the round-1 copy (or the repo source when no copy exists)
and that no two edits overlap in the same file."""
import json
import os

RUN = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z'
REPO = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
AN4 = 'RUN/analyses/AN4/20261021T124414Z'
AN6 = 'RUN/analyses/AN6/gtex-area-20261021T124444Z/out'
AN7 = 'RUN/analyses/AN7/20261021T124200Z/out'
CAT = 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out'

FAM = f'{AN4}/combined/family_summary_by_source.tsv:n_bh_q_lt_0.05 (68 GSE192741, 74 Vu); n_testable 114'
RANK = f'{AN4}/combined/focal_family_rank.tsv:rank_excess (Stromal ECM 42 GSE192741, 29 Vu)'
TOP = 'RUN/analyses/AN4/RESULT.md:18-19 (top programs by excess; hepatocyte lineage highest median excess)'
LIN = f'{AN4}/focal_lineage/focal_lineage_results.tsv rows scheme=content_lineage (excess 0.0349/0.0473, p 0.0056/0.0004 -> two-program BH q 0.0112/0.0008; ductular 0.0152 p 0.1002, -0.0011 p 0.5029)'
WCORR = f'{AN4}/focal_lineage/focal_lineage_results.tsv program_wcorr_fibroblasts vs null_wcorr_fibroblasts_mean (0.138 vs 0.090; 0.395 vs 0.263)'

edits = [
    # ---------------- spatial (AN4) ----------------
    dict(id='misc-01', file='docs/manuscript/draft/fig5.md',
         old='of that mismatch. Co-expressed programs often exceed a null of this kind, so\nexceeding it does not show program-specific tissue organization.',
         new=('of that mismatch. We therefore applied the same all-16 null, as a descriptive\n'
              'analysis that was not prespecified, to all 114 fixed programs with enough\n'
              'measured genes. In Guilliams et al., 68 of 114 programs exceeded it, and in\n'
              'Vu et al., 74 of 114 did (BH q < 0.05 over 114 programs within each source).\n'
              'ECM/IGFBP7 ranked 42nd and 29th by excess, and programs of hepatocyte\n'
              'secretory and plasma-protein genes had the largest excess. With control genes\n'
              'matched on correlation with fibroblast abundance instead, ECM/IGFBP7 still\n'
              'exceeded the null (excess 0.035, q = 0.011; 0.047, q = 0.0008), although its\n'
              'member genes remained more correlated with fibroblast abundance than the\n'
              'matched controls. With controls matched on cholangiocyte abundance,\n'
              'ductular-injury/BICC1 still did not (0.015, p = 0.10; -0.001, p = 0.50).\n'
              'Exceeding this null therefore does not show program-specific tissue\n'
              'organization.'),
         items=['B-S5'], number_sources=[FAM, RANK, TOP, LIN, WCORR], decision=None,
         reason='Replaces the round-1 qualitative "often exceed" with the measured 114-program family result and the lineage-matched null (AN4).'),
    dict(id='misc-02', file='docs/manuscript/draft/fig5.md',
         old='programs, only ECM/IGFBP7 exceeds the matched-gene null within sections, a null\nthat co-expressed programs often also exceed.',
         new=('programs, only ECM/IGFBP7 exceeds the matched-gene null within sections, but\n'
              'most fixed programs (68 and 74 of 114) also exceed it, so the spatial result\n'
              'does not show program-specific organization.'),
         items=['B-S5'], number_sources=[FAM], decision=None,
         reason='Summary sentence matches misc-01 (AN4 family result).'),
    dict(id='misc-03', file='docs/manuscript/draft/fig5_discussion.md',
         old='autocorrelated than matched control genes. This is not a donor-level result,\nand co-expressed programs often exceed the same null.',
         new=('autocorrelated than matched control genes. This is not a donor-level result,\n'
              'and most fixed programs (68 and 74 of 114 per source) exceed the same null, so\n'
              'it does not show program-specific organization.'),
         items=['B-S5'], number_sources=[FAM], decision=None,
         reason='Discussion claim matches the measured family result (AN4).'),
    dict(id='misc-04', file='docs/manuscript/draft/fig5_discussion.md',
         old='programs, and of the two only ECM/IGFBP7 exceeds a matched-gene spatial null,\nwhich co-expressed programs often also exceed.',
         new=('programs, and of the two only ECM/IGFBP7 exceeds a matched-gene spatial null,\n'
              'which most of the 114 scorable fixed programs also exceed.'),
         items=['B-S5'], number_sources=[FAM], decision=None,
         reason='Discussion claim matches the measured family result (AN4).'),
    dict(id='misc-05', file='docs/manuscript/draft/01_intro.md',
         old='a comparison that co-expressed programs often also pass.',
         new='but most of the 114 scorable fixed programs also pass this comparison, so it does not show program-specific organization.',
         items=['B-S5'], number_sources=[FAM], decision=None,
         reason='Intro claim matches the measured family result (AN4).'),
    dict(id='misc-06', file='docs/manuscript/METHODS.md',
         old='BH correction covers the two prespecified programs within each source; no prespecified test compares the two programs.',
         new=('BH correction covers the two prespecified programs within each source; no prespecified test compares the two programs. '
              'Two descriptive sensitivities were added after the primary result and were not prespecified. First, the all-16 null was applied to every fixed program with at least eight measured genes and at least 20% of its L1 weight retained (114 of 117 in each source). '
              'Each program was matched on correlation with its own discovery lineage, genes of every program were excluded from the control pool, each program had 9,999 draws with the same per-source seeds, and BH correction covered the 114 programs within each source. '
              'Second, the two focal programs were rerun with control genes matched on correlation with the abundance of their content lineage (fibroblasts for Stromal ECM, cholangiocytes for Ductular injury).'),
         items=['B-S5'], number_sources=['RUN/analyses/AN4/RESULT.md:47-50 (Task 1 method: testability rule >= 8 genes and retained L1 >= 0.20, 114 of 117, seeds 42/1042, union exclusion, BH over 114)', FAM], decision=None,
         reason='Describes the AN4 family and lineage-matched analyses that the text now reports.'),
    dict(id='misc-07', file='docs/manuscript/05_figure_legends.md',
         old='so it is not a donor-level test and is not matched for program co-expression.',
         new=('so it is not a donor-level test and is not matched for program co-expression. '
              'In a descriptive analysis that was not prespecified, 68 of 114 (Guilliams et al.) and 74 of 114 (Vu et al.) scorable fixed programs exceed the same null (BH q < 0.05 over 114 programs per source), and Stromal ECM ranks 42nd and 29th by excess. '
              'With controls matched on fibroblast abundance, Stromal ECM still exceeds the null (q = 0.011 and 0.0008); with controls matched on cholangiocyte abundance, Ductular injury does not.'),
         items=['B-S5'], number_sources=[FAM, RANK, LIN], decision=None,
         reason='Legend 5F carries the family result and lineage-matched null (AN4).'),
    dict(id='misc-08', file='docs/PAPER.md',
         old='injury does not (indeterminate); no prespecified test compared the two\n   programs.',
         new=('injury does not (indeterminate), but 68 and 74 of 114 fixed programs also\n'
              '   exceed that null, so this is not program-specific organization; no\n'
              '   prespecified test compared the two programs.'),
         items=['B-S5'], number_sources=[FAM], decision=None,
         reason='PAPER.md hero-program claim matches the AN4 family result.'),
    dict(id='misc-09', file='docs/PAPER.md',
         old="Moran's I, is therefore required before spatial organization is called\nprogram-specific.",
         new=("Moran's I, is therefore required before spatial organization is called\n"
              'program-specific. Matched-gene calibration alone is also not sufficient. In a\n'
              'descriptive analysis that was not prespecified, 68 of 114 (GSE192741) and 74\n'
              'of 114 (Vu) fixed programs exceed the same all-16-adjusted null at BH q < 0.05,\n'
              'Stromal ECM ranks 42nd and 29th by excess, and hepatocyte secretory and\n'
              'plasma-protein programs rank highest. Stromal ECM still exceeds the null with\n'
              'controls matched on fibroblast abundance, and Ductular injury still does not\n'
              'with controls matched on cholangiocyte abundance.'),
         items=['B-S5'], number_sources=[FAM, RANK, TOP, LIN], decision=None,
         reason='PAPER.md spatial finding now states the family comparison (AN4).'),
    dict(id='misc-10', file='docs/PAPER.md',
         old='Stromal ECM exceeds a matched-gene spatial null after all-cell adjustment and Ductular injury does not (indeterminate);',
         new='Stromal ECM exceeds a matched-gene spatial null after all-cell adjustment and Ductular injury does not (indeterminate), but 68 and 74 of 114 fixed programs also exceed it;',
         items=['B-S5'], number_sources=[FAM], decision=None,
         reason='Claim table row matches the AN4 family result.'),
    dict(id='misc-11', file='docs/PAPER.md',
         old='Stromal ECM exceeds its gene-level matched null in both Visium sources and Ductular injury does not;',
         new='Stromal ECM exceeds its gene-level matched null in both Visium sources and Ductular injury does not, and 68 and 74 of 114 fixed programs also exceed that null;',
         items=['B-S5'], number_sources=[FAM], decision=None,
         reason='Limitations bullet matches the AN4 family result.'),
    dict(id='misc-12', file='docs/PAPER.md',
         old='requires matched-gene calibration before it can be described\n  as program-specific organization.',
         new=('requires matched-gene calibration and comparison with the other\n'
              '  fixed programs before it can be described as program-specific\n'
              '  organization; most fixed programs exceed the matched-gene null.'),
         items=['B-S5'], number_sources=[FAM], decision=None,
         reason='Principle bullet no longer implies matched-gene calibration is sufficient (AN4).'),
    # ---------------- GTEx (AN6) ----------------
    dict(id='misc-13', file='docs/manuscript/draft/fig3.md',
         old='Only the second condition held, so the rule did not stop the analysis.',
         new=('Only the second condition held, so the rule did not stop the analysis. '
              'Both positive controls met their predictions. PLIN2, which is not a disease-state call, rose by 1.07 log2 CPM from the lowest- to the highest-fat donor (moderated *p* = 7.5 × 10⁻⁷). '
              'Across 197 steatosis-associated genes that are not disease-state calls, the GTEx fat slope correlated with the biopsy steatosis effect (rho 0.717, Freedman-Lane *p* = 0.0002). '
              'In a post hoc sensitivity, adding slide tissue area to the model left the fat correlation almost unchanged (rho 0.631 against 0.638; paired donor-bootstrap difference -0.0065, 95% interval -0.033 to 0.004).'),
         items=['PLAN-A9'], number_sources=[
             'Analysis/Hypotheses/H55-20261015T112558Z/RESULT_H55.md:25-26 (PLIN2 slope 1.069 log2 CPM, p 7.52e-07; PC3 197 genes rho 0.717, FL p 0.0002)',
             'Analysis/Hypotheses/H55-20261015T112558Z/prespec/PRESPEC_H55.md:327-336 (PLIN2 not a canonical call; PC3 = 197 H8 steatosis genes not canonical)',
             f'{AN6}/area_adjusted_summary.tsv (rho 0.631, own voom weights)',
             f'{AN6}/area_minus_primary_bootstrap.tsv (observed -0.00652, q025 -0.0327, q975 0.00409, n_boot 2000)'],
         decision=None,
         reason='Adds the AN6 area-adjusted sensitivity and the two positive controls where 0.638 is printed.'),
    dict(id='misc-14', file='docs/manuscript/METHODS.md',
         old='Negative controls were slide tissue area and ischemic time.',
         new=('Negative controls were slide tissue area and ischemic time. '
              'After the result, a post hoc sensitivity added the slide tissue-area rank as a covariate, with quality weights re-estimated on the covariate-only design that includes it. '
              'The two correlations were compared by a paired donor bootstrap (the same 2,000 resamples within the six center-by-Hardy strata, weights fixed). '
              'Resampling donors shrinks each correlation toward 0, so only the paired difference is reported.'),
         items=['PLAN-A9'], number_sources=['RUN/analyses/AN6/RESULT.md:18-21,43 (method; boot_index_Q2.rds 2,000 resamples, six strata, fixed weights; bootstrap mean 0.577 vs 0.638)'],
         decision=None,
         reason='Methods for the AN6 post hoc tissue-area sensitivity.'),
    dict(id='misc-15', file='docs/manuscript/draft/fig5_discussion.md',
         old='missed its prediction (Spearman *rho* 0.354, *p* = 0.0564) and came close to its\nstop rule.',
         new=('missed its prediction (Spearman *rho* 0.354, *p* = 0.0564) and came close to its\n'
              'stop rule. Adding tissue area to the model, a post hoc check, left the fat\n'
              'correlation at 0.631 (0.638 without it).'),
         items=['PLAN-A9'], number_sources=[f'{AN6}/area_adjusted_summary.tsv'], decision=None,
         reason='Limitation sentence now carries the area-adjusted result (AN6).'),
    dict(id='misc-16', file='docs/PAPER.md',
         old='its stop rule (rho 0.354, p 0.0564).',
         new=('its stop rule (rho 0.354, p 0.0564); adding tissue area to the model leaves\n'
              '  rho at 0.631 (post hoc; paired bootstrap difference -0.0065, 95% interval\n'
              '  -0.033 to 0.004), and PLIN2 and the 197-gene steatosis composite met their\n'
              '  positive-control predictions.'),
         items=['PLAN-A9'], number_sources=[f'{AN6}/area_adjusted_summary.tsv', f'{AN6}/area_minus_primary_bootstrap.tsv', 'Analysis/Hypotheses/H55-20261015T112558Z/RESULT_H55.md:25-26'],
         decision=None,
         reason='PAPER.md GTEx bullet matches fig3.md (AN6).'),
    # ---------------- single cell (AN7) ----------------
    dict(id='misc-17', file='docs/manuscript/draft/supplementary_notes.md',
         old='removes it (q = 0.0616 and 0.0889).',
         new=('removes it (q = 0.0616 and 0.0889).\n\n'
              '### Primary-model sensitivities\n\n'
              'Under the primary model (separate residual variance per dataset, 19,999 within-dataset stage permutations, BH over 117 programs; 55 donors for the hepatocyte programs), averaging the stored cell scores over QC-passing cells alone selects 0 of 117 programs (smallest q = 0.304). '
              'Stage betas correlate with the primary betas at Spearman 0.994. Five programs change sign, all with |beta| < 0.02 in both fits. '
              'Ductular injury gives beta 0.139 (q = 0.343) and Stromal ECM 0.076 (q = 0.578). '
              'Each stored cell score already averages over neighbors that include failing cells, so this analysis restricts the donor average and is not a rescoring on passing cells. '
              'Removing T-cell donors with fewer than 20 T cells leaves 36 of the 53 fitted T-cell donors. No T-cell program is selected (smallest q over 117 = 0.587, with the other 106 programs\' permutation p values unchanged). '
              'This is a completed test without a selection, not an untestable result.'),
         items=['B-S1', 'B-S3'], number_sources=[
             f'{AN7}/S1_summary.tsv (selected 0, smallest_perm_q 0.3042, spearman_beta_vs_h59 0.9937, sign_flips 5)',
             f'{AN7}/S1_v2qc_results_117.tsv (sign-flip programs |beta| <= 0.0168 and <= 0.0148; hepatocytes_20 beta 0.1390 q 0.3434; hepatocytes_8 beta 0.0763 q 0.5780; n_fit 55)',
             f'{AN7}/S2_summary.tsv (removed 17, n_fit 36, h59_n_fit 53, selected 0, smallest_perm_q_bh117 0.58734)'],
         decision=None,
         reason='Adds the AN7 H59-GLS QC-cell and T-cell-gate sensitivities with the not-a-rescoring caveat.'),
    dict(id='misc-18', file='docs/manuscript/METHODS.md',
         old='Projected scores are not in native Hotspot units.',
         new=('Projected scores are not in native Hotspot units. '
              'Under the primary model, two further sensitivities averaged the stored cell scores over QC-passing cells only, and removed donors with fewer than 20 T cells from the 11 T-cell program fits, with the other 106 programs\' permutation p values held fixed in the 117-program BH family (Supplementary Note 8).'),
         items=['B-S1', 'B-S3'], number_sources=['RUN/analyses/AN7/RESULT.md:100-102,136 (S2 BH over 117 with the 106 non-T-cell H59 p values unchanged)'],
         decision=None,
         reason='Methods for the AN7 sensitivities.'),
    dict(id='misc-19', file='docs/manuscript/05_figure_legends.md',
         old='17 of 53 T-cell donors have fewer than 20 T cells.',
         new='17 of 53 T-cell donors have fewer than 20 T cells; without them (36 donors), no T-cell program is selected (smallest q over 117 = 0.587).',
         items=['B-S3'], number_sources=[f'{AN7}/S2_summary.tsv (n_fit 36, selected 0, smallest_perm_q_bh117 0.58734)'],
         decision=None,
         reason='Legend 4C states the outcome of the 20-cell gate it already mentions (AN7).'),
    # ---------------- Catalog 1,853 vs 1,854 ----------------
    dict(id='misc-20', file='docs/manuscript/METHODS.md',
         old='Quarantined identities have no Catalog record and are not counted in any class.',
         new=('Quarantined identities have no Catalog record and are not counted in any class. '
              'One canonical bulk call is quarantined: ENSG00000278732 (Metazoa_SRP, a family label that eight Ensembl genes carry in the bulk fit) could not be reduced to one gene identity, so the Catalog holds 1,853 of the 1,854 calls.'),
         items=['AN5-caveat'], number_sources=[
             f'{CAT}/sidecar_ids_without_catalog_record.tsv (ENSG00000278732, GEN_IRREDUCIBLE_MULTI_ENSEMBL, bulk_canonical_call True)',
             f'{CAT}/invariants.tsv:20,26 (1854 canonical calls; canonical call without a Catalog record = ENSG00000278732)',
             'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv (8 rows with symbol Metazoa_SRP; 1,854 rows padj<0.05 & |logFC|>0.5, only ENSG00000278732 absent from catalog_r9d_gene_table.parquet rna_canonical_call)'],
         decision=None,
         reason='Explains why the Catalog (and Supp Note 7) print 1,853 while the paper prints 1,854.'),
]


def source_path(repo_file):
    rel = repo_file  # e.g. docs/manuscript/METHODS.md
    staged = os.path.join(RUN, 'staged/round1', rel)
    return staged if os.path.exists(staged) else os.path.join(REPO, rel)


problems = []
spans = {}
for e in edits:
    p = source_path(e['file'])
    s = open(p).read()
    n = s.count(e['old'])
    if n != 1:
        problems.append(f"{e['id']}: old occurs {n} times in {p}")
        continue
    i = s.index(e['old'])
    spans.setdefault(e['file'], []).append((i, i + len(e['old']), e['id']))
    for bad in ('—', 'TREAT', 'progression', 'causal'):
        if bad in e['new']:
            problems.append(f"{e['id']}: contains {bad!r}")
for f, L in spans.items():
    L.sort()
    for a, b in zip(L, L[1:]):
        if a[1] > b[0]:
            problems.append(f'overlap {a[2]} {b[2]} in {f}')

notes = (
    'Round-2 writer misc. All old strings checked unique against RUN/staged/round1 copies (build script analyses/R2/misc/build_main_misc.py). '
    'No jobs were submitted (no compute needed; job range 3231-3239 unused). '
    'Catalog spatial sidecar to stage: Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN4/20261021T124414Z/combined/catalog_spatial_rows_all16_sidecar.tsv '
    '(4 Visium rows keyed by program_context_id; all-16 estimates, P, two-program q, 114-program rank/percentile/q, content-lineage excess/q; no state change from r9d). '
    'Abstract has no spatial claim, so no abstract edit. '
    'Ancestry Holm (item 3): round-1 bulk-20 (fig3.md) and bulk-21 (Supp Note 7) are correct against AN6 (one of 200 / one of 100 reference fits at or above observed; 2/201 = 0.00995, 2/101 = 0.0198; Holm 0.0199) and neither says "floor". No ancestry "floor" wording remains in any round-1 copy; only the bulk-20/21 reason fields (not paper text) still say "resolution floor". No edit. '
    'The 1,853 vs 1,854 gap: ENSG00000278732 (Metazoa_SRP; logFC 0.535, padj 1.6e-5) is a canonical call quarantined as GEN_IRREDUCIBLE_MULTI_ENSEMBL in Catalog r9d, so it has no Catalog record (r9d invariants.tsv rows 20, 26, 51). Supp Note 7 (cat-22) already says the Catalog holds 1,853 of the calls; misc-20 adds the reason in Methods. The same gene is the one 233-sidecar gene without an r9d record (AN6). '
    'misc-01 uses "Guilliams et al." for GSE192741, as fig5.md does. Family values for the two focal programs differ slightly from the two-program primary because of the larger exclusion set (AN4 Task 1); the text prints only family counts and ranks from the family run, and lineage-matched values from the focal run. '
    'Item IDs: PLAN-A9 is the GTEx area adjustment (critic MISSED item in C7); AN5-caveat is the 1,853/1,854 trace.'
)

out = dict(edits=edits, panels=[], notes=notes, problems=problems)
with open(os.path.join(RUN, 'staged/main_misc.json'), 'w') as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print('edits', len(edits), 'problems', problems)
