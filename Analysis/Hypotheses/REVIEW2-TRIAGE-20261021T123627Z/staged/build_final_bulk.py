import json, copy, re
REPO='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
RUN=REPO+'Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/'
src=json.load(open(RUN+'staged/edits_bulk.json'))['edits']
E={e['id']:copy.deepcopy(e) for e in src}
AN3='Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/'
def mod(i, **kw):
    E[i].update(kw)

# bulk-01: clarify exclusion rule; split last sentence
mod('bulk-01', new=("Genotypes were called from the RNA-seq reads at 8,660 autosomal sites, keeping calls with read depth of at least 10. "
 "Libraries with fewer than 1,000 called sites were not used. Two libraries were joined into one genetic-identity cluster when KING-robust kinship (PLINK 2.0) was at least 0.354 over at least 500 shared sites and also in each of two disjoint halves of the sites (at least 250 shared sites in each half). "
 "A cluster was excluded when any candidate same-person pair in it failed this rule or when its libraries carried conflicting source labels. "
 "Kinship was computed for every pair of the 1,466 libraries from ten cohorts in one joint panel, so the screen covered pairs within and across studies. "
 "All 72 pairs at or above 0.354 were within one cohort. Each pair of the 774 selected libraries shared at least 598 sites, and no pair reached 0.354. "
 "One pair from different cohorts, both disease units (GSE213621 and GSE130970), had kinship 0.288 (0.294 and 0.282 in the two halves) and opposite homozygous genotypes at 54 of 2,845 shared sites, against at most 6 in any pair at or above 0.354. "
 "This fits first-degree relatives, not one person, and both libraries were kept. No other pair of selected libraries reached kinship 0.0884 in PLINK. "
 "Because no pair from different cohorts reached 0.354, the screen has no positive example of a cross-cohort duplicate. KING also cannot distinguish monozygotic twins from repeated libraries of one person."),
 number_sources=E['bulk-01']['number_sources']+["checker re-derivation 2026-10-21: bcftools query of identity_panel.merged.vcf.gz at the 8,660 filtered_identity_site_selection.tsv sites, DP>=10, 774 bulk_primary_sample_ids.tsv: min pairwise shared sites 598 over 299,151 pairs; min called sites per library 1,003",
 "full_filtered_related_pairs.tsv: 114 rows (PLINK --king-table-filter 0.0884); only SRR21622926/SRR9036348 has both libraries in the 774 (re-derived)"])

# bulk-04: clarify sign wording
mod('bulk-04', new=E['bulk-04']['new'].replace("Of the 1,854 primary calls, 1,779 are also called there, and none has the opposite sign.",
  "Of the 1,854 primary calls, 1,779 are also called there, all with the same sign."))
assert "all with the same sign" in E['bulk-04']['new']

# bulk-06: D3 user answer, five cohorts with GSE126848 filterByExpr refit
mod('bulk-06', old=("We next fit each cohort separately with the same 24,115-gene family and canonical threshold (Figure 3C). The union of per-cohort calls contains 13,761 genes: 9,065 called in one cohort and 4,696 in two or more. Only 2,831 of the 4,696 have the same sign in every cohort that calls them (1,933 in two cohorts, 636 in three, 218 in four and 44 in all five). The other 1,865 (39.7%) are called up in one cohort and down in another, and 1,632 of these include a GSE126848 call. Without GSE126848, 346 of 2,565 multi-cohort genes (13.5%) have mixed signs. These 20,054 cohort memberships are a different set from the 1,854 pooled calls, and shared membership is not replication."),
 new=("We next fit each cohort separately with the canonical threshold (Figure 3C). Four cohorts use the 24,115-gene family of the pooled fit. "
 "Because GSE126848 calls track sequencing depth (see below), we refit that cohort on the 15,592 genes that pass its own expression filter (filterByExpr), with BH correction over those 15,592 genes. The refit gives 3,440 calls. "
 "The union of per-cohort calls contains 12,180 genes: 8,341 called in one cohort and 3,839 in two or more. Only 2,628 of the 3,839 have the same sign in every cohort that calls them (1,820 in two cohorts, 582 in three, 186 in four and 40 in all five). "
 "The other 1,211 (31.5%) are called up in one cohort and down in another, and 946 of these include a GSE126848 call. Without GSE126848, 346 of 2,565 multi-cohort genes (13.5%) have mixed signs. "
 "With the 6,147 GSE126848 calls from the 24,115-gene family in place of the refit, 1,865 of 4,696 multi-cohort genes (39.7%) have mixed signs. "
 "These 17,347 cohort memberships are a different set from the 1,854 pooled calls, and shared membership is not replication."),
 number_sources=[AN3+"fig3c_union_filterByExpr_refit.tsv rows 'GSE126848 within-cohort filterByExpr refit' (17347, 12180, 8341, 3839, 1211, 31.54, 946, same2-5 1820/582/186/40), 'four cohorts without GSE126848' (2565, 346, 13.49) and 'adopted calls (five cohorts)' (4696, 1865, 39.71)",
  AN3+"gse126848_filterByExpr_accounting.tsv (filterByExpr_genes 15592, filterByExpr_calls 3440)",
  AN3+"an3_gse126848_filtered_union.R:31-33 (filterByExpr, RLE, voomWithQualityWeights on filtered genes; topTable adj.P.Val = BH over 15,592)",
  "checker re-derivation 2026-10-21 from figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv + gse126848_filterByExpr_refit_genes.tsv.gz: all union counts reproduced exactly; BH over 15,592 reproduced (max diff 1.4e-15); 2,628 = 1,820+582+186+40"],
 decision="D3",
 reason="Implements the user's D3 answer (DECISIONS.md): five cohorts with GSE126848 entering through its within-cohort filterByExpr refit; states the 15,592 vs 24,115 BH-family rule. Replaces the four-cohort form.")

# bulk-08: legend in D3 form
mod('bulk-08', old=("**(C)** Per-cohort calls (BH FDR < 0.05 and |log2FC| > 0.5 within each cohort's 24,115-gene family). Top, calls by direction. Middle, genes per overlap class in the 13,761-gene union; darker bars include a GSE126848 call. Bottom, mixed-sign share of multi-cohort genes, with and without GSE126848. Shared membership is not replication. GSE126848 calls (7 disease, 23 control clusters) also track library depth, and GSE135251"),
 new=("**(C)** Per-cohort calls (BH FDR < 0.05 and |log2FC| > 0.5). Four cohorts use the fixed 24,115-gene family. GSE126848 calls come from a refit on the 15,592 genes that pass its own expression filter, with BH over those 15,592 genes. "
 "Top, calls by direction. Middle, genes per overlap class in the 12,180-gene union; darker bars include a GSE126848 call. Bottom, mixed-sign share of multi-cohort genes: 1,211 of 3,839 (31.5%) with GSE126848 and 346 of 2,565 (13.5%) without it. "
 "Shared membership is not replication. In GSE126848 (7 disease, 23 control clusters), disease libraries have about half the depth of controls, and GSE135251"),
 number_sources=[AN3+"fig3c_union_filterByExpr_refit.tsv", AN3+"gse126848_filterByExpr_accounting.tsv (15592)", "docs/manuscript/draft/fig3.md:9 (median 19.5 versus 35.3 million fragments)"],
 decision="D3",
 reason="Fig 3C legend in the user's D3 form; must match a five-cohort filtered-refit redraw, which does not exist yet (the AN3 candidate is four-cohort).")

# bulk-09 / bulk-10 PAPER.md D3 form
mod('bulk-09', new=("GSE126848 enters per-cohort membership through a refit on the 15,592 genes that pass its own expression filter (BH over those genes; the other four cohorts use the fixed 24,115-gene family). "
 "The per-cohort union contains 12,180 genes, including 8,341 called in one cohort and 3,839 in two or more, of which 2,628 have one sign in every calling cohort (40 in all five) and 1,211 (31.5%) have opposite signs in different cohorts; it is distinct from the 1,854 pooled calls."),
 number_sources=[AN3+"fig3c_union_filterByExpr_refit.tsv row 'GSE126848 within-cohort filterByExpr refit'", AN3+"gse126848_filterByExpr_accounting.tsv (15592)"],
 decision="D3", reason="PAPER.md Fig 3 row in the user's D3 form (five cohorts, GSE126848 filtered refit).")
mod('bulk-10', new="shared per-cohort membership is not replication, because 1,211 of 3,839 multi-cohort genes change sign between cohorts.",
 number_sources=[AN3+"fig3c_union_filterByExpr_refit.tsv row 'GSE126848 within-cohort filterByExpr refit' (3839, 1211)"],
 decision="D3", reason="PAPER.md claim column in the user's D3 form.")

# bulk-11: same-sign wording for GSE162694
mod('bulk-11', new=E['bulk-11']['new'].replace("and the next cohort, GSE162694, calls 1,014.", "and the next highest, GSE162694, calls 1,014 with the same sign."),
 number_sources=E['bulk-11']['number_sources']+[AN3+"pooled_1854_called_by_each_cohort_alone.tsv (GSE213621 1654 same sign; GSE162694 1015 overlap, 1014 same sign)"])
assert "next highest" in E['bulk-11']['new']
E['bulk-13']['number_sources']=E['bulk-13']['number_sources']+[AN3+"pooled_1854_called_by_each_cohort_alone.tsv"]

# bulk-14 heading: drop the F0-reference clause
mod('bulk-14', new="## Expression differs from F0 at F2 to F4, and the F1 difference depends mainly on one cohort",
 reason="'Every fibrosis stage' overstated F1 (14 of 1,338 without GSE130970). F0-reference clause dropped: with a disease-only F0, F1 still gives 783 genes (bulk-15), so 'depends on' overstates; no library-property clause per D4.")

# bulk-15: wording
mod('bulk-15', new=E['bulk-15']['new'].replace("F1 gives 783 genes", "F1 versus F0 gives 783 genes"))

# bulk-20 / 21: tie to D5; tighten wording
mod('bulk-20', new=E['bulk-20']['new'].replace("so both p values are limited by the number of reference fits.", "so the precision of both p values is limited by the number of reference fits."),
 decision="D5", reason=E['bulk-20']['reason']+" Holds only if D5 keeps the B = 200 / 100 nulls (no B = 1,000 rerun).")
mod('bulk-21', decision="D5", reason=E['bulk-21']['reason']+" Holds only if D5 keeps the B = 200 / 100 nulls.")

# bulk-23: wording
mod('bulk-23', new=("Memberships, weights and coordinates were built on an earlier 23,370-gene expression universe and are unchanged. "
 "System scores use the 24,115 genes tested in the 774-unit data. Because 392 genes of the earlier universe are absent from these data, four nodes keep less than 80% of their members and are untestable."))
# bulk-24: avoid duplicating bulk-25's 'map was built' sentence
mod('bulk-24', new=("The system redundancy summary used all 466 GSE162694 and GSE213621 units (140 and 326), including their 31 and 68 controls, and is descriptive. "
 "Because 392 genes of the 23,370-gene universe on which the map was built are absent from the 24,115-gene family, four nodes keep less than 80% of their original membership after signature exclusion and are untestable."))
# bulk-26: wording
mod('bulk-26', new=("Of the 145 signature genes, 143 are observed in the 24,115-gene family; DPEP1 and ZACN are not. "
 "In the earlier 23,370-gene family used for the molecular-systems map and the expression-relationship extension, 139 were observed, and KRT23, PRAMEF10, IQGAP3 and IER3 were also absent. The mapping table lists the unavailable genes."))

# bulk-29 + bulk-30 merged
mod('bulk-29', old="(k6 NMF; the sixth axis, k6-P6, is unavailable). These axes describe continuous variation, not patient classes.",
 new=("(k6 NMF). These axes describe continuous variation, not patient classes. The sixth axis, k6-P6, is not shown. "
 "With its fixed weights restricted to the 4,517 retained genes, its scores recomputed from the original counts correlated with the released scores at Spearman 0.9874, below the prespecified minimum of 0.99. "
 "On the fragment counts used here, the five shown axes correlate with the released scores at Spearman 0.928 (k6-P5) to 0.998 (k6-P1) over 617 shared libraries."),
 number_sources=E['bulk-29']['number_sources']+E['bulk-30']['number_sources'],
 reason="Gives the k6-P6 reason (historical-count reconstruction, per C8 critic) and the corrected-count agreement of the five shown axes; merged with bulk-30 so the two Spearman values are not read as the same check.")
del E['bulk-30']

# bulk-31: avoid two adjacent parentheticals
mod('bulk-31', old="Bottom, F0-centered activity of five k6 NMF axes (continuous, not patient classes)",
 new="Bottom, F0-centered activity of five k6 NMF axes (continuous, not patient classes; k6-P6 fails the fixed reconstruction rule and is not shown)")

# bulk-32: version wording
mod('bulk-32', new=("The signal-level direction analysis used an earlier version of the 50-study colocalization results (2026-08-17; 462 SuSiE genes), "
 "not the corrected rerun that gives the 436 SuSiE genes reported elsewhere in this paper."))

# new D3 edits
E['bulk-33']={"id":"bulk-33","file":"docs/manuscript/METHODS.md",
 "old":("Figure 3C shows the unadjusted `~ group_binary` calls. Three sensitivities are reported in the text (Supplementary Note 7) and are not used for the panel: refits of `~ group_binary + log2(library size)` with library size counted before or after RLE normalization; the number of calls on genes with zero reads in all 30 libraries; and a refit after `filterByExpr` within the 30 units with RLE factors recomputed (3,440 calls on 15,592 genes, a different BH family)."),
 "new":("Figure 3C uses a refit after `filterByExpr` within the 30 units, with RLE factors recomputed and BH over the 15,592 retained genes, a different family from the fixed 24,115 genes of the other four cohorts (3,440 calls). "
 "Two sensitivities of the 6,147 calls are reported in the text (Supplementary Note 7) and are not used for the panel: refits of `~ group_binary + log2(library size)` with library size counted before or after RLE normalization, and the number of calls on genes with zero reads in all 30 libraries."),
 "items":["B-B2"],"number_sources":[AN3+"gse126848_filterByExpr_accounting.tsv (15592, 3440)"],"decision":"D3",
 "reason":"METHODS.md:133 says Figure 3C shows the unfiltered calls and the filterByExpr refit is not used for the panel; under the user's D3 answer the panel uses the refit."}
E['bulk-34']={"id":"bulk-34","file":"docs/manuscript/draft/supplementary_notes.md",
 "old":"A cohort-specific expression filter gives 3,440 calls on 15,592 genes.",
 "new":"A cohort-specific expression filter gives 3,440 calls on 15,592 genes, with BH over those genes, and Figure 3C uses these calls.",
 "items":["B-B2"],"number_sources":[],"decision":"D3",
 "reason":"Supplementary Note 7 names the GSE126848 calls that Figure 3C uses under the user's D3 answer."}

order=[e['id'] for e in src if e['id'] in E]+['bulk-33','bulk-34']
final=[E[i] for i in order]

# checks
from collections import defaultdict
spans=defaultdict(list); errs=[]
banned=['—','TREAT','causal','progression','sealed provenance','promotion','release gate','validator','audit trail','staging']
for e in final:
    t=open(REPO+e['file']).read()
    n=t.count(e['old'])
    if n!=1: errs.append((e['id'],'count',n))
    i=t.find(e['old']); spans[e['file']].append((i,i+len(e['old']),e['id']))
    for b in banned:
        if b in e['new'] and b not in e['old']: errs.append((e['id'],'banned',b))
for f,s in spans.items():
    s.sort()
    for a,b in zip(s,s[1:]):
        if b[0]<a[1]: errs.append(('overlap',f,a[2],b[2]))
print('errors',errs, 'n',len(final))
json.dump({"final_edits":final}, open('/scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/39a23967-7721-4ff8-86d5-d2d5894e47f5/scratchpad/final_edits_only.json','w'), indent=1)
# render applied copies for reading
out={}
for f,s in spans.items():
    t=open(REPO+f).read()
    for e in sorted([x for x in final if x['file']==f], key=lambda x: t.find(x['old']), reverse=True):
        i=t.find(e['old']); t=t[:i]+e['new']+t[i+len(e['old']):]
    out[f]=t
import pickle; pickle.dump(out, open('/scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/39a23967-7721-4ff8-86d5-d2d5894e47f5/scratchpad/applied.pkl','wb'))

# ---- second pass ----
E['bulk-15'].update(old="but only GSE130970 gives many calls on its own.",
 new="but only GSE130970 gives many calls on its own. With F0 restricted to its 89 disease biopsies (191 units), F1 versus F0 gives 783 genes (Figure S3E).",
 reason="Reports the disease-only F0 count where the F1 dependence is discussed; placed after the cohort sentences so the cohort argument is not split.")
E['bulk-18']['new']=E['bulk-18']['new'].replace(
 "This section is a post hoc descriptive deviation from the prespecified reading of the library-property check, whose positive control failed (Methods).",
 "This section qualifies the fibrosis and NAS stage counts (Figures 3E and S3D). It is a post hoc descriptive deviation from the prespecified reading of the library-property check, whose positive control failed (Methods).")
assert "qualifies the fibrosis and NAS stage counts" in E['bulk-18']['new']
final=[E[i] for i in order]
spans=defaultdict(list); errs=[]
for e in final:
    t=open(REPO+e['file']).read()
    n=t.count(e['old'])
    if n!=1: errs.append((e['id'],'count',n))
    i=t.find(e['old']); spans[e['file']].append((i,i+len(e['old']),e['id']))
    for b in banned:
        if b in e['new'] and b not in e['old']: errs.append((e['id'],'banned',b))
for f,s in spans.items():
    s.sort()
    for a,b in zip(s,s[1:]):
        if b[0]<a[1]: errs.append(('overlap',f,a[2],b[2]))
# cross-writer overlap
import glob
for fpath in glob.glob(RUN+'staged/*.json'):
    if 'bulk' in fpath: continue
    d=json.load(open(fpath)); rows=d if isinstance(d,list) else d.get('final_edits',d.get('edits',[]))
    for o in rows:
        try: t=open(REPO+o['file']).read()
        except Exception: continue
        i=t.find(o['old'])
        for e in final:
            if e['file']!=o['file'] or i<0: continue
            j=t.find(e['old'])
            if i<j+len(e['old']) and j<i+len(o['old']): errs.append(('xoverlap',fpath.split('/')[-1],o['id'],e['id']))
print('second pass errors',errs,'n',len(final))
FINAL=final

dropped=[{"id":"bulk-30","why":"Merged into bulk-29 (same fig3.md:25 passage) so the 0.9874 reconstruction rule and the 0.928-0.998 corrected-count agreement are stated together and are not read as the same check."}]
added=["bulk-33","bulk-34"]
problems=[
 {"id":"bulk-06/08/09/10","kind":"decision_conflict","evidence":"Writer staged the four-cohort Fig 3C form (10,630; 13.5%). DECISIONS.md D3 chose five cohorts with GSE126848 via its filterByExpr refit. Rewritten from Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/fig3c_union_filterByExpr_refit.tsv (12,180; 3,839 multi; 1,211 mixed, 31.5%; 946 with GSE126848; same-sign 1,820/582/186/40; 17,347 memberships). Checker re-derived every count from cohort_disease_all_gene_results.tsv + gse126848_filterByExpr_refit_genes.tsv.gz, and BH over 15,592 genes (max diff 1.4e-15). Tags kept as D3 but these edits implement the user's answer; do not drop them with --drop-decision D3."},
 {"id":"bulk-33","kind":"contradiction_missed","evidence":"METHODS.md:133 said 'Figure 3C shows the unadjusted ~ group_binary calls' and listed the filterByExpr refit as 'not used for the panel'. Under D3 that is false. Added bulk-33 (METHODS.md:133) and bulk-34 (supplementary_notes.md:73)."},
 {"id":"fig3c-candidate","kind":"figure_mismatch","evidence":"The only Fig 3C candidate render, figures/candidates/fig3c-four-cohort-20261021T124839Z (AN3 jobs 3122-3126, caption '...four cohorts without GSE126848 (10,630 genes)'), is in the four-cohort form that D3 rejected. Legend bulk-08 describes a five-cohort filtered-refit panel (12,180-gene union; 31.5% and 13.5% bottom shares; GSE126848 top bars should show 3,440 calls, 1,821 up / 1,619 down). A D3 redraw is needed before install; read its printed values against bulk-08."},
 {"id":"bulk-14","kind":"overstatement","evidence":"Heading 'the F1 difference depends on one cohort and on the F0 reference': with a disease-only F0, F1 still gives 783 genes (H82 d4_f0_reference.tsv base_diseaseF0: 783 calls, 678 of 1,338 retained), so 'depends on' overstates. Heading now drops 'every fibrosis stage' and keeps only the cohort clause; bulk-15 reports 783 in the body."},
 {"id":"bulk-31","kind":"readability","evidence":"Writer's legend insertion produced two adjacent parentheticals ('five k6 NMF axes (k6-P6 ...) (continuous, not patient classes)'). Merged into one parenthetical."},
 {"id":"bulk-20/21","kind":"hold_rule","evidence":"Ancestry p-value text depends on PLAN A8 / D5 (rerun at B = 1,000 is the alternative). Tagged D5 so they drop if the null is rerun. Values verified: RESULT_H79.md:21-22 (#ref >= obs 1 and 1; p 0.00995, 0.0198; Holm 0.0199)."},
 {"id":"bulk-01","kind":"wording","evidence":"'A cluster that contained an unresolved pair' was undefined; rewritten as 'any candidate same-person pair in it failed this rule' per 11_resource_identity.py:128-152, 308-315. The 598 minimum shared sites (299,151 pairs) was re-derived by the checker with bcftools at DP >= 10, not taken from the critic."},
 {"id":"plan-A3.1","kind":"input_error","evidence":"PLAN.md section 2 gives the cross-cohort pair 'IBS0 27'. 27 is the per-half value (half0/half1_filtered_related_pairs.tsv); the full-site value is 54 of 2,845 (full_filtered_related_pairs.tsv), which bulk-01 prints."},
 {"id":"C8-critic-31","kind":"input_error","evidence":"C8 critic calls the 31 GSE162694 units 'unstaged'. model_design.tsv (bulk-corrected-controls-20261008T162308Z) gives 31 controls in GSE162694 and 68 in GSE213621; bulk-22/24 say controls, which is correct."},
 {"id":"D4-install","kind":"outside_targets","evidence":"bulk-16/17/18 report H82 counts as a post hoc deviation. RESULT_H82.md and docs/RESULTS.md:1472-1475 say no H82 number enters the paper; a DEVIATIONS.md entry is needed (draft at analyses/AN3/h82/DEVIATIONS_entry_draft.md) and RESULTS.md must change at install."}
]
summary=("33 final edits (32 staged minus bulk-30 merged, plus bulk-33/34). All old strings occur once; no overlaps within bulk or with edits_cat/gen/sc. "
 "Main change: D3 edits rewritten from the four-cohort form to the user's five-cohort filterByExpr-refit form (12,180 union; 1,211 of 3,839 mixed, 31.5%), with the 15,592 vs 24,115 BH-family rule stated in fig3.md, legend, METHODS and PAPER; METHODS.md:133 contradiction fixed. "
 "Fig 3E heading no longer claims F0-reference dependence. Re-derived: identity 72/0 cross-cohort, 598 min shared sites, pair 0.288/54 of 2,845; Fig 3C unions; LOCO 0.997-1.000; 1,654/1,014; H82 11-contrast counts; H56 controls; H79 p values; 466 = 140+326 with 31+68 controls; k6 Spearman values; strict-depth 738/24,509/1,906/1,779. "
 "Fig 3C still needs a D3 redraw before the legend can be installed.")
json.dump({"final_edits":FINAL,"dropped":dropped,"added_ids":added,"problems":problems,"summary":summary},
          open(RUN+'staged/final_bulk.json','w'), indent=1, ensure_ascii=False)
print('saved', RUN+'staged/final_bulk.json')
