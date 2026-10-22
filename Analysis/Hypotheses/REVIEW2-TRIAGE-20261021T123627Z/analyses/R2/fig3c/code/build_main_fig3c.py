"""Round-2 fig3c writer: build RUN/staged/main_fig3c.json and check every 'old' occurs once in the round-1 copy
(or in the repo source when no copy exists), and that no two edits overlap in one file."""
import json, os
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
RUN = 'Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/'
R1 = ROOT + RUN + 'staged/round1/'
AN3 = RUN + 'analyses/AN3/run-20261021T124631Z/'
ACCT = AN3 + 'gse126848_filterByExpr_accounting.tsv'
UNION = AN3 + 'fig3c_union_filterByExpr_refit.tsv'
CAND = 'figures/candidates/fig3c-five-filtered-20261021T172748Z/'

E = []
def add(i, f, old, new, items, srcs, dec, reason):
    E.append(dict(id=i, file=f, old=old, new=new, items=items, number_sources=srcs, decision=dec, reason=reason))

add('r2-fig3c-01', 'docs/manuscript/draft/fig3.md',
    'Its model contains 23 control and 7 disease clusters, all male, and gives 6,147 calls.',
    'Its model on the fixed 24,115-gene family contains 23 control and 7 disease clusters, all male, and gives 6,147 calls.',
    ['B-B2'], [], 'D3',
    'Paragraph 7 now gives 3,440 refit calls for Figure 3C; this names which fit gives the 6,147 so the two counts are not read as conflicting.')

add('r2-fig3c-02', 'docs/manuscript/draft/supplementary_notes.md',
    'A cohort-specific expression filter gives 3,440 calls on 15,592 genes, with BH over those genes, and Figure 3C uses these calls.',
    'A refit after the cohort\'s own expression filter (`filterByExpr` within the 30 units, RLE factors recomputed) gives 3,440 calls (1,821 up, 1,619 down) on 15,592 genes, with BH over those genes, and Figure 3C uses these calls. All 3,258 refit calls that are also among the 6,147 keep their sign, and 182 refit calls are new.',
    ['B-B2'], [ACCT + ' (filterByExpr_up 1821, filterByExpr_down 1619, refit_calls_shared_with_adopted 3258, refit_calls_shared_same_sign 3258, refit_calls_not_in_adopted 182)',
               CAND + 'check/fig3c_cohort_overlap_bars.txt line 2 (3,440 (1,821 up, 1,619 down))'], 'D3',
    'The panel prints the GSE126848 refit split 1,821 up / 1,619 down; the note now gives the same split and how the refit calls relate to the 6,147.')

add('r2-fig3c-03', 'docs/PAPER.md',
    'and none of the 437, so the GSE126848 calls are not read as disease replication.',
    'and none of the 437, so the GSE126848 calls are not read as disease replication.\nFigure 3C uses a refit on the 15,592 genes that pass the cohort\'s own expression\nfilter (BH over those genes; 3,440 calls, 1,821 up and 1,619 down).',
    ['B-B2'], [ACCT + ' (filterByExpr_genes 15592, filterByExpr_calls 3440, up 1821, down 1619)'], 'D3',
    'The PAPER.md GSE126848 paragraph described only the 6,147-call fit; under D3 the panel uses the filtered refit.')

add('r2-fig3c-04', 'docs/manuscript/05_figure_legends.md',
    'darker bars include a GSE126848 call.',
    'darker segments are genes with a GSE126848 call.',
    ['B-B2'], [], 'D3',
    'In the drawn panel each class bar is split into two segments; the dark segment, not the whole bar, counts genes with a GSE126848 call.')

add('r2-fig3c-05', 'docs/STATUS.md',
    '- Decision 8: Figure 3C shows the mixed-sign share of multi-cohort genes with and without GSE126848 (39.7% and 13.5%).',
    '- Decision 8 (revised 2026-10-21, D3): Figure 3C uses five cohorts, with GSE126848 entering through its within-cohort filterByExpr refit (BH over its 15,592 genes; the other four cohorts use the fixed 24,115-gene family). It shows the mixed-sign share of multi-cohort genes with and without GSE126848 (31.5% and 13.5%).',
    ['B-B2'], [UNION + " rows 'GSE126848 within-cohort filterByExpr refit' (31.54) and 'four cohorts without GSE126848' (13.49)", ACCT + ' (15592)'], 'D3',
    'Apply at install of the five-filtered 3C panel: Decision 8 fixed 39.7%, which the user replaced under D3.')

add('r2-fig3c-06', 'docs/STATUS.md',
    'The five-cohort union contains 13,761 genes: 9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts; 1,865 of the 4,696 multi-cohort genes have opposite signs between cohorts (1,632 involve GSE126848), leaving 1,933 / 636 / 218 / 44 single-sign genes in two to five cohorts.',
    'Since 2026-10-21 (decision D3), Figure 3C uses the GSE126848 within-cohort filterByExpr refit instead (3,440 calls; BH over its 15,592 genes). The five-cohort union contains 12,180 genes: 8,341 called in one cohort and 3,839 in two or more; 1,211 of the 3,839 have opposite signs between cohorts (946 involve GSE126848), leaving 1,820 / 582 / 186 / 40 single-sign genes in two to five cohorts. With the 6,147 calls the union was 13,761 genes and 1,865 of 4,696 (39.7%) were mixed.',
    ['B-B2'], [UNION + " row 'GSE126848 within-cohort filterByExpr refit' (12180, 8341, 3839, 1211, 946, same2-5 1820/582/186/40); row 'adopted calls (five cohorts)' (13761, 4696, 1865, 39.71)", ACCT + ' (15592, 3440)'], 'D3',
    'Apply at install: the H1 row described the 13,761-gene union as the current Figure 3C.')

add('r2-fig3c-07', 'docs/STATUS.md',
    'Its 6,147 calls join a 13,761-gene per-cohort union (9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts; 1,865 of the 4,696 multi-cohort genes have opposite signs between cohorts), distinct from pooled membership.',
    'Figure 3C uses its within-cohort filterByExpr refit (3,440 calls, BH over 15,592 genes; decision D3, 2026-10-21), which joins a 12,180-gene per-cohort union (8,341 called in one cohort and 3,839 in two or more; 1,211 of the 3,839 have opposite signs between cohorts), distinct from pooled membership.',
    ['B-B2'], [UNION + " row 'GSE126848 within-cohort filterByExpr refit'", ACCT + ' (15592, 3440)'], 'D3',
    'Apply at install: the bulk row described the 6,147-call union as the current Figure 3C.')

add('r2-fig3c-08', 'docs/ROADMAP.md',
    'endpoint roster 392, five-cohort union 13,761',
    'endpoint roster 392, five-cohort union 12,180 (GSE126848 filtered refit, D3 2026-10-21)',
    ['B-B2'], [UNION + " row 'GSE126848 within-cohort filterByExpr refit' (12180)"], 'D3',
    'Apply at install: checklist row printed the superseded 3C union.')

add('r2-fig3c-09', 'docs/RESULTS.md',
    'to one cohort) and the historical 844-record censuses are superseded for this panel.',
    'to one cohort) and the historical 844-record censuses are superseded for this panel.\n\n'
    '**Figure 3C with the GSE126848 filtered refit, 2026-10-21 (user decision D3).**\n'
    'GSE126848 enters per-cohort membership through a refit after `filterByExpr` within\n'
    'its 30 units (RLE factors recomputed; BH over its 15,592 genes; 3,440 calls, 1,821\n'
    'up and 1,619 down). The other four cohorts keep the fixed 24,115-gene family. The\n'
    'union contains 12,180 genes and 17,347 memberships: 8,341 genes called in one\n'
    'cohort and 3,839 in two or more. 2,628 have one sign in every calling cohort\n'
    '(1,820 / 582 / 186 / 40 in two to five cohorts) and 1,211 (31.5%) are mixed, 946 of\n'
    'them with a GSE126848 call. Without GSE126848, 346 of 2,565 (13.5%) are mixed; with\n'
    'the 6,147 calls, 1,865 of 4,696 (39.7%). All 3,258 refit calls that are also among\n'
    'the 6,147 keep their sign, and 182 refit calls are new. 2,730 of the 6,147 calls\n'
    '(2,639 up) are on genes the cohort\'s own filter removes. Sources:\n'
    '`Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/`\n'
    '(job 3121, seed 20260923) and `figures/candidates/fig3c-five-filtered-20261021T172748Z/`.\n'
    'This replaces the 2026-10-14 union above for Figure 3C.',
    ['B-B2'], [UNION + " rows 'GSE126848 within-cohort filterByExpr refit' (17347, 12180, 8341, 3839, 1211, 31.54, 946, 1820/582/186/40), 'four cohorts without GSE126848' (2565, 346, 13.49), 'adopted calls (five cohorts)' (4696, 1865, 39.71)",
               ACCT + ' (15592, 3440, 1821, 1619, 3258, 182, adopted_calls_on_genes_removed_by_filter 2730, adopted_calls_removed_by_filter_up 2639)',
               '2,628 = 1,820 + 582 + 186 + 40'], 'D3',
    'Apply at install: RESULTS.md records the new 3C numbers with source paths, appended after the 2026-10-14 block, which stays as history.')

add('r2-fig3c-10', 'docs/RESULTS.md',
    'reading no H82 number enters the paper, and the review numbers on library complexity of\nthe stage contrasts stay review-only. Figure 3E counts are unchanged.',
    'reading no H82 number enters the paper. On 2026-10-21 the user approved a logged post hoc\ndeviation (decision D4; `H82-20261015T125052Z/DEVIATIONS.md`, "Post hoc reporting\ndeviation"): the D1 joint-adjustment counts for all 11 primary-test contrasts are reported\ndescriptively in Supplementary Note 7, with no F1-specific wording and no T1/T2 result\n(table: `REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/h82/h82_joint_adjustment_retention_11_contrasts.tsv`).\nFigure 3E counts are unchanged.',
    ['A3.2', 'B-B4'], [], 'D4',
    'Apply with the H82 DEVIATIONS entry: RESULTS.md said no H82 number enters the paper, which the D4 supplementary line contradicts.')

out = []
for e in E:
    p = R1 + e['file'] if os.path.exists(R1 + e['file']) else ROOT + e['file']
    t = open(p).read()
    n = t.count(e['old'])
    assert n == 1, (e['id'], n, p)
    e['_start'] = t.index(e['old']); e['_end'] = e['_start'] + len(e['old']); e['_path'] = p
for a in E:
    for b in E:
        if a is not b and a['file'] == b['file']:
            assert a['_end'] <= b['_start'] or b['_end'] <= a['_start'], (a['id'], b['id'])
for e in E:
    print(e['id'], e['file'], 'checked in', 'round1 copy' if e['_path'].startswith(R1) else 'repo source')
    for k in ['_start', '_end', '_path']: del e[k]
    assert '—' not in e['new'], e['id']
panels = [dict(candidate_pdf=CAND + 'panels/fig3c_cohort_overlap_bars.pdf',
               install_target='figures/main/fig3_bulk_transcriptomics/panels/fig3c_cohort_alluvial.pdf')]
json.dump(dict(edits=E, panels=panels), open(ROOT + RUN + 'staged/main_fig3c.json', 'w'), indent=1, ensure_ascii=False)
print(len(E), 'edits written')
