"""Blind checker for round-2 writer misc. Reads the writer's staged/main_misc.json,
applies checker corrections, verifies each old string occurs exactly once in the
round-1 copy (or repo source), checks no overlap with any other staged/main_*.json,
simulates applying the edits to scratch copies, and overwrites staged/main_misc.json
as {"final_edits": [...], "panels": []}."""
import glob
import json
import os
import sys

RUN = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z'
REPO = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
OUT = RUN + '/staged/main_misc.json'
SIM = sys.argv[1] if len(sys.argv) > 1 else None

AN4 = 'RUN/analyses/AN4/20261021T124414Z'
AN6 = 'RUN/analyses/AN6/gtex-area-20261021T124444Z/out'
CAT = ('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/'
       'catalog_r9d-20261016T105550Z/build_out')
DEG = ('RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/'
       'primary_source_controls/deg_results.csv')

WRITER = '/scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/39a23967-7721-4ff8-86d5-d2d5894e47f5/scratchpad/main_misc.writer.json'
src = json.load(open(WRITER))
writer = {e['id']: e for e in (src.get('final_edits') or src.get('edits'))}

TOP10 = (f'{AN4}/combined/family_all_programs.tsv: top 10 by excess per source (GSE192741: C3 acute-phase, '
         'APOA1 plasma protein x2, APOA2, SAA1, ALB, APOE, APOC1, FOS, F12; Vu: CYP xenobiotic, ALB, C3, APOA2, '
         'GOT1, SAA1, APOA1 x2, CFH, CRP)')
LIN_Q = (f'{AN4}/focal_lineage/focal_lineage_results.tsv scheme=content_lineage residual_padj_two_program '
         '(Stromal ECM 0.0112/0.0008; Ductular injury 0.1002/0.5029); excess 0.0349/0.0473 and 0.0152/-0.0011')

fix = {}

fix['misc-01'] = dict(
    new=('of that mismatch. We therefore applied the matched-gene null, under the\n'
         'all-16 model and with each program\'s controls matched on its own discovery\n'
         'lineage, to all 114 fixed programs with enough measured genes, as a\n'
         'descriptive analysis that was not prespecified. In Guilliams et al., 68 of 114\n'
         'programs exceeded it, and in Vu et al., 74 of 114 did (BH q < 0.05 over 114\n'
         'programs within each source). ECM/IGFBP7 ranked 42nd and 29th of 114 by\n'
         'excess, and programs of hepatocyte plasma-protein, acute-phase and metabolic\n'
         'genes had the largest excess. With control genes matched on correlation with\n'
         'fibroblast abundance instead, ECM/IGFBP7 still exceeded the null (excess\n'
         '0.035, q = 0.011; 0.047, q = 0.0008), although its member genes remained more\n'
         'correlated with fibroblast abundance than the matched controls. With controls\n'
         'matched on correlation with cholangiocyte abundance, ductular-injury/BICC1\n'
         'still did not (0.015, q = 0.10; -0.001, q = 0.50). Exceeding this null\n'
         'therefore does not show program-specific tissue organization.'),
    extra_sources=[TOP10, LIN_Q],
    reason_suffix=(' Checker: stated that each program is matched on its own discovery lineage; "matched on '
                   'cholangiocyte abundance" corrected to "on correlation with"; Vu top program is xenobiotic CYP, '
                   'so "secretory and plasma-protein" widened to "plasma-protein, acute-phase and metabolic"; '
                   'Ductular given as two-program q like ECM.'))

fix['misc-05'] = dict(
    new=('but most of the 114 scorable fixed programs also pass this comparison, so the result does not show '
         'program-specific organization.'),
    reason_suffix=' Checker: "so it does not" had an unclear referent.')

fix['misc-06'] = dict(
    new=('BH correction covers the two prespecified programs within each source; no prespecified test compares '
         'the two programs. Two descriptive analyses with the 16-abundance model were added after the primary '
         'result and were not prespecified. First, the matched-gene null was applied to every fixed program with '
         'at least eight measured genes and at least 20% of its L1 weight retained (114 of 117 in each source). '
         'Control genes for each program were matched on correlation with the abundance of its discovery lineage, '
         'genes of every fixed program were excluded from the control pool, each program had 9,999 draws from one '
         'random stream per source (seeds 42 and 1042), and BH correction covered the 114 programs within each '
         'source. Second, the two prespecified programs were rerun with control genes matched on correlation with '
         'the abundance of the lineage in which each scores highest (fibroblasts for Stromal ECM, cholangiocytes '
         'for Ductular injury), with BH over the two programs.'),
    extra_sources=['RUN/analyses/AN4/RESULT.md:47-50 (testability rule, 114 of 117, own-lineage matching, '
                   'one RNG per source seeds 42/1042, union exclusion, BH over 114); RESULT.md:22-31 (content-lineage '
                   'run, q is BH over two programs)'],
    reason_suffix=(' Checker: removed the undefined terms "all-16 null" and "content lineage"; added that the '
                   'lineage-matched rerun uses BH over two programs.'))

fix['misc-07'] = dict(
    new=('so it is not a donor-level test and is not matched for program co-expression. In a descriptive analysis '
         'that was not prespecified, 68 of 114 (Guilliams et al.) and 74 of 114 (Vu et al.) scorable fixed programs '
         'exceed a null built the same way, with each program\'s controls matched on its own discovery lineage '
         '(BH q < 0.05 over 114 programs per source); Stromal ECM ranks 42nd and 29th by excess. With controls matched on '
         'correlation with fibroblast abundance, Stromal ECM still exceeds the null (q = 0.011 and 0.0008); with '
         'controls matched on correlation with cholangiocyte abundance, Ductular injury does not (q = 0.10 and 0.50).'),
    extra_sources=[LIN_Q],
    reason_suffix=(' Checker: the family null matches each program on its own lineage, not on hepatocyte abundance, '
                   'so "the same null" was inexact; "matched on fibroblast abundance" corrected to "on correlation '
                   'with"; Ductular q added.'))

fix['misc-09'] = dict(
    new=('Moran\'s I, is therefore required before spatial organization is called\n'
         'program-specific. Matched-gene calibration alone is also not sufficient. In a\n'
         'descriptive analysis that was not prespecified, 68 of 114 (GSE192741) and 74\n'
         'of 114 (Vu) fixed programs exceed a matched-gene null built the same way, with\n'
         'all-16 adjustment and controls matched on each program\'s own discovery\n'
         'lineage (BH q < 0.05 over 114 programs per source). Stromal ECM ranks 42nd and\n'
         '29th by excess, and programs of hepatocyte plasma-protein, acute-phase and\n'
         'metabolic genes rank highest. Stromal ECM still exceeds the null with controls\n'
         'matched on correlation with fibroblast abundance (q=0.011; q=0.0008), and\n'
         'Ductular injury still does not with controls matched on correlation with\n'
         'cholangiocyte abundance (q=0.10; q=0.50).'),
    extra_sources=[TOP10, LIN_Q],
    reason_suffix=(' Checker: same corrections as misc-01 (own-lineage matching, "correlation with", top-program '
                   'description, q values).'))

fix['misc-13'] = dict(
    new=('Only the second condition held, so the rule did not stop the analysis. Both positive controls met their '
         'predictions. PLIN2, which is not a disease-state call, rose by 1.07 log2 CPM from the lowest- to the '
         'highest-fat donor (moderated *p* = 7.5 × 10⁻⁷). Across 197 steatosis-associated genes that are not '
         'disease-state calls, the GTEx fat slope correlated with the biopsy steatosis effect (rho 0.717, '
         'Freedman-Lane *p* = 0.0002, the smallest two-sided value 10,000 permutations allow). In a post hoc '
         'sensitivity, adding slide tissue area to the model left the fat correlation almost unchanged (rho 0.631 '
         'against 0.638; difference -0.0065, paired donor-bootstrap 95% interval -0.033 to 0.004).'),
    extra_sources=[f'{AN6}/area_adjusted_summary.tsv primary_reproduced pc3_fl_p_two 0.00019998 = 2/10,001 '
                   '(two-sided floor; RESULT_H55.md:32 "floor 2/10001" for the same draws)'],
    reason_suffix=(' Checker: PC3 Freedman-Lane p is two-sided and at the 2/10,001 floor (AN6 column pc3_fl_p_two), '
                   'so the floor is stated; the bootstrap gives the interval, not the -0.0065 point difference.'))

fix['misc-14'] = dict(
    new=('Negative controls were slide tissue area and ischemic time. After the result, a post hoc sensitivity '
         'added the slide tissue-area rank as a covariate, with quality weights re-estimated on the covariate-only '
         'design that includes it. The area-adjusted and primary correlations were compared on the same 2,000 '
         'stratified donor-bootstrap resamples (six center-by-Hardy strata, weights fixed). Donor resampling '
         'shifts each correlation toward 0 (bootstrap mean 0.577 against the observed 0.638), so only the paired '
         'difference is reported.'),
    extra_sources=[f'{AN6}/bootstrap_rho.tsv.gz (mean rho_primary 0.5768 over 2,000 resamples)'],
    reason_suffix=' Checker: "shrinks" replaced by the observed bootstrap shift with its number.')

fix['misc-15'] = dict(
    new=('missed its prediction (Spearman *rho* 0.354, *p* = 0.0564) and came close to its\n'
         'stop rule. In a post hoc check, adding tissue area to the model left the fat\n'
         'correlation at 0.631 (0.638 without it).'),
    reason_suffix=' Checker: sentence order simplified.')

fix['misc-17'] = dict(
    new=writer['misc-17']['new'].replace('averages over neighbors that include failing cells',
                                         'averages over neighbors that can include failing cells'),
    reason_suffix=' Checker: "include" softened to "can include" (neighbor sets vary by cell; AN7 RESULT.md:58-60).')

fix['misc-18'] = dict(
    new=('Projected scores are not in native Hotspot units. Two further sensitivities used the primary model. One '
         'averaged the stored cell scores over QC-passing cells only. The other removed donors with fewer than 20 '
         'T cells from the 11 T-cell program fits and kept the other 106 programs\' permutation p values in the '
         '117-program BH family (Supplementary Note 8).'),
    reason_suffix=' Checker: split one run-on sentence.')

fix['misc-20'] = dict(
    old='The Catalog therefore holds 1,853 of the 1,854 disease-versus-control RNA calls.',
    new=('One disease-versus-control call, ENSG00000278732, is quarantined because its symbol, Metazoa_SRP, is a '
         'family label that eight Ensembl IDs carry in the bulk fit. The Catalog therefore holds 1,853 of the '
         '1,854 disease-versus-control RNA calls.'),
    sources=[f'{CAT}/sidecar_ids_without_catalog_record.tsv:86 (ENSG00000278732, GEN_IRREDUCIBLE_MULTI_ENSEMBL, '
             'Metazoa_SRP, bulk_canonical_call True)',
             f'{CAT}/invariants.tsv (bulk canonical calls 1854; canonical call without a Catalog record = '
             'ENSG00000278732)',
             f'{DEG} (24,115 rows; 1,854 with padj < 0.05 and |logFC| > 0.5; 8 distinct gene IDs with symbol '
             'Metazoa_SRP, 1 of them a call)'],
    reason_suffix=(' Checker: moved to the sentence that already states 1,853 of 1,854 (METHODS quarantine '
                   'paragraph), which names the multi-Ensembl quarantine rule; the writer\'s placement repeated '
                   'the 1,853 count two sentences later.'))

fix['misc-16'] = dict(
    new=writer['misc-16']['new'].replace('(post hoc; paired bootstrap difference -0.0065, 95% interval',
                                         '(post hoc; difference -0.0065, paired bootstrap 95% interval'),
    reason_suffix=' Checker: the bootstrap gives the interval, not the point difference (as misc-13).')

final = []
for eid, e in writer.items():
    e = dict(e)
    f = fix.get(eid)
    if f:
        if 'old' in f:
            e['old'] = f['old']
        e['new'] = f['new']
        if 'sources' in f:
            e['number_sources'] = f['sources']
        e['number_sources'] = list(e.get('number_sources', [])) + f.get('extra_sources', [])
        e['reason'] = e['reason'] + f.get('reason_suffix', '')
    final.append({k: e.get(k) for k in ['id', 'file', 'old', 'new', 'items', 'number_sources', 'decision', 'reason']})

# checks: exactly once, no overlap within misc and against other main_*.json
def text_of(f):
    p = RUN + '/staged/round1/' + f
    return open(p if os.path.exists(p) else REPO + '/' + f).read()

spans = {}
ok = True
for e in final:
    t = text_of(e['file'])
    c = t.count(e['old'])
    if c != 1:
        ok = False
        print('COUNT', e['id'], c)
    i = t.find(e['old'])
    spans.setdefault(e['file'], []).append((i, i + len(e['old']), 'misc', e['id']))
    for bad in ['—', 'TREAT', 'causal', 'progression']:
        if bad in e['new'] and bad not in e['old']:
            print('LANG', e['id'], bad)
for jf in glob.glob(RUN + '/staged/main_*.json'):
    if jf.endswith('main_misc.json'):
        continue
    d = json.load(open(jf))
    for e in (d.get('final_edits') or d.get('edits') or []):
        t = text_of(e['file'])
        i = t.find(e['old'])
        spans.setdefault(e['file'], []).append((i, i + len(e['old']), os.path.basename(jf), e['id']))
for f, s in spans.items():
    s.sort()
    for a, b in zip(s, s[1:]):
        if b[0] < a[1]:
            ok = False
            print('OVERLAP', f, a, b)
print('checks ok' if ok else 'CHECKS FAILED', len(final))

if SIM:
    by = {}
    for e in final:
        by.setdefault(e['file'], []).append(e)
    for f, es in by.items():
        t = text_of(f)
        for e in es:
            t = t.replace(e['old'], e['new'], 1)
        os.makedirs(os.path.dirname(os.path.join(SIM, f)), exist_ok=True)
        open(os.path.join(SIM, f), 'w').write(t)

if ok:
    json.dump({'final_edits': final, 'panels': []}, open(OUT, 'w'), indent=1, ensure_ascii=False)
    print('wrote', OUT)
