# Build staged/final_gen.json from staged/edits_gen.json with checker corrections.
import json, os, copy

ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
RUN = os.path.join(ROOT, 'Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z')
src = json.load(open(os.path.join(RUN, 'staged/edits_gen.json')))
edits = {e['id']: copy.deepcopy(e) for e in src['edits']}
order = [e['id'] for e in src['edits']]

def fix(eid, reason_add, **kw):
    e = edits[eid]
    e.update(kw)
    e['reason'] = e['reason'] + ' Checker: ' + reason_add

# gen-intro-1
fix('gen-intro-1',
    'Seagle 2025 abstract says genetically predicted gene expression, not liver expression; PMID placeholders removed because final_cat cat-77 adds Hong 2025 (ref 42) and Seagle 2025 (ref 45), and Zhu 2026 is ref 40.',
    new=("Multi-ancestry GWAS have been combined with genetically predicted gene expression to prioritize MASLD genes (Seagle et al., 2025). "
         "Single-cell eQTL mapping in liver biopsies from 25 people with MASLD and 23 controls found cell-state-dependent eQTLs enriched for MASLD heritability (Hong et al., 2025). "
         "Liver chromatin accessibility, reporter-assay, liver eQTL and CRISPR interference data have been combined to assign target genes to MASLD risk variants (Zhu et al., 2026). "
         "None of these studies keeps colocalization of fine-mapped multi-ancestry GWAS with a European liver-eQTL panel, multicohort disease-state expression and donor-aware cell and tissue context side by side for each gene."))

# gen-fig2-2
fix('gen-fig2-2',
    'Removed the AST/ALT/NAFLD split because two genes (WHRN, MAU2) tie within MVP EUR and the split depends on tie order; the total 140 does not. Reworded the internal-fit sentence, which implied one shared fit.',
    new=("Of the same 436 genes, 94 are supported only by the Million Veteran Program (MVP) European ALT, AST, and NAFLD GWAS [ref: MVP liver-trait GWAS], which had no standalone fine-mapping, "
         "and these three GWAS supply the strongest colocalization for 140 of the 436 genes. "
         "SuSiE was refit to the GWAS inside each gene-by-GWAS colocalization test, in the same way for every stratum, whether or not the stratum was fine-mapped on its own (Methods)."))

# gen-fig2-5
fix('gen-fig2-5',
    'TOR1B is not tied: MVP_NAFLD_EUR PP.H4 0.9999927 against PDFF 0.9999840, so the parenthesis is removed.',
    new=("12 direct only, 11 of them SuSiE-tested in an enzyme GWAS, and 10 with both classes. "
         "MVP European NAFLD is the only direct-trait support for 9 of the 22 genes and gives the strongest direct-trait PP.H4 for 11 of the 22."))

# gen-fig2-9 (old extended to the end of the sentence)
fix('gen-fig2-9',
    'The stated reason (wrong positions in identifiers) is false: norm_key strips trailing zeros, so truncation cannot cause a miss. Per the C5 critic, 44 of the 82 UKBB credible sets were fit at an earlier lead set. Added why the colocalization-targeted UKBB run (loci selected at PP.H4 >= 0.5, combined_finemapping_coloc_targeted.csv) is excluded. Old extended to the sentence end so the sentence closes cleanly.',
    old=("These GWAS were fine-mapped at 6, 9, and 15 loci, which match 8 of their 767 lead windows by locus identifier (23 of 767 within 500 kb by position), "
         "and the arm has no credible set at PNPLA3 or TM6SF2."),
    new=("In the fine-mapping table used for Figure 2D, these GWAS have credible sets at 6, 9, and 15 loci and none at PNPLA3 or TM6SF2. "
         "These loci came from an earlier lead set, so they match 8 of the 767 current lead windows by locus identifier, while 23 of 767 have a credible set within 500 kb by position, the better coverage measure. "
         "A separate UKBB fine-mapping run at loci selected by colocalization (PP.H4 ≥ 0.5) is not used, because loci selected on colocalization cannot serve as an unbiased enzyme arm."),
    number_sources=edits['gen-fig2-9'].get('number_sources', []) + [
        'inputs/C5-finemap.json:285 (44 of 82 UKBB credible sets not at a current lead; earlier lead set)',
        'inputs/C5-finemap.json:202-203,287 (targeted run selected at COLOC PP.H4 >= 0.5; 00c_prep_coloc_targeted_leads.R:5-13)'])

# gen-supp-1
fix('gen-supp-1',
    'Added that every alternative interval includes 1 (RESULT_H57.md:39-51); "twofold enrichment" changed to "an odds ratio of 2"; "rests on one matched list" stated as its cause (all 1,000 matching orders returned the same list).',
    new=("The deciles did not fully control expression, so we read no direction from this result. "
         "The upper limit of the interval, 2.007, does not exclude an odds ratio of 2. "
         "Other adjustments gave odds ratios from 0.958 to 1.347, all with 95% intervals that include 1: 5 and 20 expression bins 1.297 and 1.304, a spline in expression 1.315 (95% CI 0.847 to 2.042), expression measured in control units only 1.347 (95% CI 0.867 to 2.093), standard-error deciles 1.028, expression by standard-error bins 0.958, and one gene per 500-kb cluster 1.288. "
         "The matched expression check (13 against 21.3 expected) used a single matched list, because all 1,000 matching orders returned the same list. "
         "Expression curves fitted to non-colocalized genes predict 20.2 to 20.8 disease-state calls among the 426 colocalized genes, close to the 21.2 expected."))

# gen-supp-2 (absorbs dropped cat-12 opening, per final_cat)
fix('gen-supp-2',
    'Denominator corrected per C1 critic: 10,935 = 10,928 without a multi-signal result + 7 SuSiE-untestable (fig2.md:21 prints 10,928). Added "95% CI" to the ABF interval. Opening sentence replaced with the cat-12 wording that final_cat dropped in favour of this edit (Catalog call indeterminate, PLAN B-C2). D5 and D8 do not depend on PLAN A1: A1 can only turn carried negatives untestable, and D5 already counts those as non-colocalized; D8 uses ABF.',
    new=("These genes have no multi-signal colocalization result. Their Catalog genetic call is indeterminate, and they are not genetic negatives. "
         "Over all 14,704 shared genes, counting the 10,935 genes outside the evaluated set (10,928 without a multi-signal result and 7 SuSiE-untestable) as non-colocalized, 26 genes were in both sets against 22.3 expected within expression deciles (odds ratio 1.19, 95% CI 0.79 to 1.81). "
         "With the single-signal ABF map in place of SuSiE (PP.H4 ≥ 0.5), 49 genes were in both against 43.6 expected (odds ratio 1.15, 95% CI 0.85 to 1.57)."),
    items=['A1.1', 'B-C2'],
    number_sources=edits['gen-supp-2'].get('number_sources', []) + [
        'inputs/C1-overlap.json:111-112 (14,704 = no_susie_result 10,928 + tested_le_0.5 3,343 + positive 426 + untestable 7)'])

# gen-supp-4
fix('gen-supp-4',
    'The European control failed only at the 98.333% level (N1a 181 of 200 at 95% against a floor of 180; 188 against 190 at 98.333%; N1_coverage.tsv, fig2.md:39). No new counts printed, because PLAN A2 will recompute these controls.',
    new=("The European, AMR and combined estimates all failed their random-split calibration controls. "
         "In these controls the intervals excluded zero more often than a calibrated interval would. "
         "The European interval failed only at the Bonferroni-corrected 98.333% level, but the AMR and combined intervals failed already at 95%, so these two intervals are too narrow. "
         "The pooled direct-minus-enzyme difference comes mainly"))

# gen-supp-5
fix('gen-supp-5',
    'Restoration rule stated as implemented (C5 B-G4.2): the lead position is read from the locus identifier and dropped trailing zeros are restored by the one restoration within 1 Mb of all member variants.',
    new=("To group credible sets into regions, we placed each set at the lead position in its locus identifier and joined sets by single linkage when their lead positions lay within 1 Mb, so one region can span more than 1 Mb. "
         "Where the identifier had lost trailing zeros, we restored the one position within 1 Mb of all the set's member variants."))

# gen-meth-2 (old extended to remove "fully")
fix('gen-meth-2',
    'Old extended to the sentence end so "fully ancestry-matched" becomes "ancestry-matched", which an all-European panel makes redundant. Does not overlap gen-meth-1.',
    old=("The Broadaway liver-eQTL panel is predominantly European [ref: Broadaway et al., 2024], so we describe analyses that pair non-European GWAS with this panel as cross-population rather than fully ancestry-matched colocalization."),
    new=("The Broadaway liver-eQTL panel is entirely of European ancestry [ref: Broadaway et al., 2024], so we describe analyses that pair non-European GWAS with this panel as cross-population rather than ancestry-matched colocalization."))

# gen-meth-4
fix('gen-meth-4',
    'n is N_eff only for binary traits and N_tot for quantitative traits (06_susie_coloc.R:265-270; eligibility-run logs print [binary:N_eff] and [quant:N_tot]). SuSiEx ran at all 81 loci and returned at least one credible set at 77 (4 .cs files are NULL); MESuSiE returned a fit at 79. SuSiEx --pval_thresh filters credible sets with no marginal p below 1e-5 (SuSiEx --help); it is not a locus threshold. MESuSiE coverage 0.95 is the meSuSie_get_cs default. "Neither fit sets a random seed" removed (deterministic fits).',
    new=("The GWAS-side SuSiE fit inside the colocalization analysis was refit for each gene-by-GWAS pair on the GWAS variants at the positions of the gene's eQTL fit, with no p-value threshold, so it ran for every stratum whether or not that stratum had standalone fine-mapping. "
         "It required at least 100 GWAS-eQTL overlapping variants and at least 50 variants with LD at the eQTL-fit positions, used the ancestry-matched GWAS LD matrix with 1 × 10⁻³ added to the diagonal (1 × 10⁻⁶ in standalone fine-mapping), z scores with n equal to the effective sample size for binary traits and the total sample size for quantitative traits, L = 10, a minimum absolute within-set correlation of 0.5, at most 500 iterations and residual variance fixed at 1, and the seed was set to 42 before each fit.\n\n"
         "Cross-ancestry fine-mapping within MVP used SuSiEx [version to be supplied] [ref: SuSiEx] and MESuSiE 1.0 [ref: MESuSiE]. "
         "Shared loci were defined per MVP trait where a lead with p < 5 × 10⁻⁸ in at least two MVP ancestries (EUR, AFR, AMR, EAS) lay within 500 kb, and the window spanned the contributing leads ± 100 kb. "
         "For ALT, AST and NAFLD this gave 81 loci. MESuSiE returned a fit at 79 of them, and SuSiEx returned at least one credible set at 77. "
         "SuSiEx computed LD from 1000 Genomes genotypes for all four ancestries and used credible-set coverage 0.95, minimum purity 0.5, MAF ≥ 0.005 and its multi-step choice of the number of effects, and it removed credible sets with no variant at marginal p < 1 × 10⁻⁵. "
         "MESuSiE used PolyFun UK Biobank LD for EUR and 1000 Genomes LD for the other ancestries with 1 × 10⁻³ added to the diagonal, L = 10, at most 200 iterations, residual variance not estimated, minimum absolute within-set correlation 0.5 and coverage 0.95, and arms with fewer than 50 variants were dropped."),
    number_sources=[
        'GWAS/finemapping/src/06_susie_coloc.R:111-120,265-270,712-735',
        'GWAS/finemapping/logs/coloc_eligibility/coloc_elig_*.out ([binary:N_eff], [quant:N_tot])',
        'GWAS/finemapping/src/09b_identify_shared_loci_mvp.R:45-48,160-205',
        'GWAS/finemapping/src/10_run_susiex.py:441-447; bin/SuSiEx --help (--pval_thresh filters credible sets)',
        'GWAS/finemapping/src/10b_run_mesusie.R:36-38,380-382; MESuSiE::meSuSie_get_cs coverage = 0.95',
        'GWAS/finemapping/src/run_susiex_mvp.sh, run_mesusie_mvp.sh (LD panels)',
        'GWAS/finemapping/results/susiex_mvp/shared_loci.csv (ALT 43, AST 35, NAFLD 3 = 81)',
        'GWAS/finemapping/results/mesusie_mvp/mesusie_locus_summary_mvp.csv (ALT 42, AST 34, NAFLD 3 = 79)',
        'GWAS/finemapping/results/susiex_mvp/{ALT,AST,NAFLD}/*.cs (81 run, 4 NULL: 77 with a credible set)',
        'micromamba/envs/rnaseq/lib/R/library/MESuSiE/DESCRIPTION Version 1.0'])

# gen-meth-5
fix('gen-meth-5',
    'Number source corrected: the r2 sessionInfo.txt is the assembler session and does not list coloc; coloc 5.2.3 is in susie_coloc_eligibility_20260923T211525Z/env/batch_1_sessionInfo.txt and the first-run logs. Priors at 06_susie_coloc.R:111-113.',
    number_sources=[
        'GWAS/finemapping/results/susie_coloc_eligibility_20260923T211525Z/env/batch_1_sessionInfo.txt (coloc_5.2.3, susieR_0.14.2)',
        'GWAS/finemapping/logs/coloc_* ("This is coloc version 5.2.3")',
        'GWAS/finemapping/src/06_susie_coloc.R:111-113 (p1 = p2 = 1e-4, p12 = 5e-6)'])

# gen-meth-8
fix('gen-meth-8',
    'Per C4: 94 of 108 genes have at least one own-refit row, so "the GWAS-caQTL entry" becomes "every qualifying GWAS-caQTL entry" for 14. NARS2: the largest -log10 p in its credible set is 0.74, so the smallest p is 0.18; MTMR3 is a one-variant set with p 0.096.',
    new=("A region was linked to a Gene Catalog gene only when, at the same peak, the caQTL row also colocalized (coloc.bf_bf PP.H4 > 0.8, same priors) with the eQTL SuSiE effect of that gene that had colocalized with the same GWAS lead variant. "
         "The eQTL fits store GRCh37 positions without alleles, so these positions were lifted to GRCh38 with the same chain file and matched to caQTL positions on position only; where a position held more than one caQTL variant, the variant with the smallest p was matched. "
         "The lead-variant LD rule was not applied to this comparison, so the 108 linked genes are an upper count. "
         "The GWAS signal was matched across GWAS refits by its lead variant, so for 14 of the 108 genes every qualifying GWAS-caQTL entry comes from another gene's refit that shares the same GWAS lead. "
         "For two genes, NARS2 and MTMR3, the colocalized eQTL SuSiE effect has no marginal eQTL association in the Broadaway summary statistics: the smallest marginal p in the credible set is 0.18 for NARS2, and the one variant in the MTMR3 set has p = 0.096."))

# gen-leg-3
fix('gen-leg-3',
    'UKBB AST does have a credible set at EFHD1 (chr2:233,509,761, susie_cs 1) in the colocalization-targeted run (combined_finemapping_coloc_targeted.csv); credible_sets.csv has no UKBB_AST chr2 locus. Both the original and the first draft were wrong as absolute statements.',
    new='UKBB AST has no credible set at this locus in the untargeted fine-mapping run, although the colocalization-targeted run has one')

added = [
    dict(id='gen-fig2-16', file='docs/manuscript/draft/fig2.md',
         old='The MVP European ALT, AST, and NAFLD GWAS were not fine-mapped.',
         new='The MVP European ALT, AST, and NAFLD GWAS had no standalone fine-mapping.',
         items=['B-G4.4'], number_sources=[], decision=None,
         reason='Checker addition: matches gen-fig2-2 and gen-meth-4, which state that colocalization refit SuSiE to these GWAS; "were not fine-mapped" would contradict them.'),
    dict(id='gen-leg-4', file='docs/manuscript/05_figure_legends.md',
         old='the MVP European ALT, AST and NAFLD strata were not fine-mapped',
         new='the MVP European ALT, AST and NAFLD strata had no standalone fine-mapping',
         items=['B-G4.4'], number_sources=[], decision=None,
         reason='Checker addition: same contradiction as gen-fig2-16, in the Figure S2A legend.'),
    dict(id='gen-paper-5', file='docs/PAPER.md',
         old='pooled 1-Mb clustering yields 233 clusters',
         new='pooled single-linkage clustering of lead positions within 1 Mb yields 233 clusters',
         items=['B-G4.2'], number_sources=[], decision=None,
         reason='Checker addition: PAPER.md:615 still says "1-Mb clustering"; clusters are chained by single linkage (largest span 6.18 Mb). No number changed.'),
]
for a in added:
    edits[a['id']] = a
    order.append(a['id'])

final = [edits[i] for i in order]

# ---- validation: once-only and no overlap within gen and against other writers' final edits
cache = {}
def text(f):
    if f not in cache:
        cache[f] = open(os.path.join(ROOT, f)).read()
    return cache[f]
errs = []
spans = {}
for e in final:
    s = text(e['file'])
    n = s.count(e['old'])
    if n != 1:
        errs.append(f"{e['id']}: old count {n}")
        continue
    i = s.find(e['old'])
    spans.setdefault(e['file'], []).append((i, i + len(e['old']), 'gen', e['id']))
for w in ('bulk', 'cat', 'sc'):
    p = os.path.join(RUN, f'staged/final_{w}.json')
    for e in json.load(open(p))['final_edits']:
        s = text(e['file'])
        if s.count(e['old']) != 1:
            continue
        i = s.find(e['old'])
        spans.setdefault(e['file'], []).append((i, i + len(e['old']), w, e['id']))
for f, ss in spans.items():
    ss.sort()
    for a in range(len(ss)):
        for b in range(a + 1, len(ss)):
            if ss[b][0] < ss[a][1] and ('gen' in (ss[a][2], ss[b][2])):
                errs.append(f'overlap {f}: {ss[a]} {ss[b]}')
for e in final:
    if '—' in e['new'] or '–' in e['new'].replace('3112–3125', ''):
        pass
    if '—' in e['new']:
        errs.append(f"{e['id']}: em dash")
print('edits', len(final), 'errors', errs)

out = dict(
    final_edits=final,
    dropped=[],
    added_ids=[a['id'] for a in added],
)
json.dump(out, open(os.path.join(RUN, 'scratch/check_gen/final_gen_edits_only.json'), 'w'), indent=1, ensure_ascii=False)

# ---- second pass: 8-of-767 wording elsewhere
more = [
    dict(id='gen-fig2-17', file='docs/manuscript/draft/fig2.md',
         old='coverage (enzyme arm fine-mapped at 8 of 767 lead windows, direct arm at 47 of 47)',
         new='coverage (enzyme credible sets at 8 of 767 lead windows by locus identifier and 23 by position, direct arm at 47 of 47)',
         items=['B-G4.3'], number_sources=['Analysis/Hypotheses/H58-20261015T112421Z G7_coverage.tsv (EUR enzyme 767 / 8 / 23), per inputs/C5-finemap.json:37,231'], decision=None,
         reason='Checker addition: after gen-fig2-9, "fine-mapped at 8 of 767" would contradict the 30 fine-mapped loci and the statement that 23 of 767 is the better measure. EUR leads are clumped at 5e-8, so PLAN A2 does not change these counts.'),
    dict(id='gen-supp-7', file='docs/manuscript/draft/supplementary_notes.md',
         old='(the UKBB enzyme GWAS were fine-mapped at 8 of 767 lead windows)',
         new='(UKBB enzyme credible sets lie at 8 of 767 lead windows by locus identifier and within 500 kb of 23)',
         items=['B-G4.3'], number_sources=['G7_coverage.tsv (EUR enzyme 767 / 8 / 23), per inputs/C5-finemap.json:37,231'], decision=None,
         reason='Checker addition: same contradiction as gen-fig2-17, in Supplementary Note 4.'),
]
for a in more:
    edits[a['id']] = a
    order.append(a['id'])
final = [edits[i] for i in order]
errs = []
spans = {}
for e in final:
    s = text(e['file'])
    n = s.count(e['old'])
    if n != 1:
        errs.append(f"{e['id']}: old count {n}"); continue
    i = s.find(e['old'])
    spans.setdefault(e['file'], []).append((i, i + len(e['old']), 'gen', e['id']))
for w in ('bulk', 'cat', 'sc'):
    for e in json.load(open(os.path.join(RUN, f'staged/final_{w}.json')))['final_edits']:
        s = text(e['file'])
        if s.count(e['old']) != 1: continue
        i = s.find(e['old'])
        spans.setdefault(e['file'], []).append((i, i + len(e['old']), w, e['id']))
for f, ss in spans.items():
    ss.sort()
    for a in range(len(ss)):
        for b in range(a + 1, len(ss)):
            if ss[b][0] < ss[a][1] and ('gen' in (ss[a][2], ss[b][2])):
                errs.append(f'overlap {f}: {ss[a]} {ss[b]}')
for e in final:
    if '—' in e['new']: errs.append(e['id'] + ' em dash')
print('pass2 edits', len(final), 'errors', errs)
out = dict(final_edits=final, dropped=[], added_ids=[a['id'] for a in added + more])
json.dump(out, open(os.path.join(RUN, 'scratch/check_gen/final_gen_edits_only.json'), 'w'), indent=1, ensure_ascii=False)
