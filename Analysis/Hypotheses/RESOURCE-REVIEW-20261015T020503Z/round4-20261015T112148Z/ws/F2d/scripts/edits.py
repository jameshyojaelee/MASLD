import json
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
E=[]
def add(fid,file,old,new,claim=False):
    E.append(dict(fid=fid,file=file,old_text=old,new_text=new,claim_level=claim))
AB='docs/manuscript/draft/abstract.md'; IN='docs/manuscript/draft/01_intro.md'; F2='docs/manuscript/draft/fig2.md'
LG='docs/manuscript/05_figure_legends.md'; DI='docs/manuscript/draft/fig5_discussion.md'; ME='docs/manuscript/METHODS.md'
PA='docs/PAPER.md'; ST='docs/STATUS.md'; RE='docs/RESULTS.md'
# ---------- R13-08 framing
add('R13-08',AB,"""Inherited genetic variation and the molecular state of established disease
capture different dimensions of metabolic dysfunction-associated steatotic
liver disease (MASLD), but""","""Inherited genetic variation and the molecular state of established disease
have each been measured in metabolic dysfunction-associated steatotic
liver disease (MASLD), but""",True)
add('R13-08',F2,"## Trait-associated inherited variation defines a genetic map distinct from established disease-state transcription",
"## Trait-associated inherited variation, qualified by phenotype, ancestry and fine-mapping coverage",True)
add('R13-08',F2,"""The distinction between inherited genetic evidence and established tissue
state has two implications. First,""","""The colocalized and disease-state gene sets overlap at about the rate expected
for sets of these sizes (26 shared against 33.46 expected); this is a
description, not evidence that the two layers are biologically separate. Two
practical points follow. First,""")
add('R13-08',F2,"""colocalized liver-trait regulatory candidates. It complements, but does not
substitute for, the established disease-state transcriptome.""","""colocalized liver-trait regulatory candidates. It does not substitute for the
established disease-state transcriptome, and the two gene sets overlap at about
the rate expected for sets of these sizes.""")
# ---------- R02-01 / R14-02 / R14-06 coding ratio: abstract
add('R02-01',AB,"""Credible sets for
direct MASLD and liver-fat phenotypes carried
23.34% of their posterior probability in protein-altering or canonical-splice
variants, against 3.87% for the liver enzymes ALT, AST and GGT, a 6.026671-fold
ratio (95% physical-cluster bootstrap interval 2.574593 to 13.510165;
difference 0.194695, 0.069229 to 0.339940). Across 57 leave-one-study-out,
leave-one-cluster-out, window-size and LD-block comparisons the ratio ranged
from 4.368132 to 8.275691 and every
difference interval stayed above zero.""","""Pooled across
ancestries, credible sets for direct MASLD and liver-fat phenotypes carried
23.34% of their posterior probability in protein-altering or canonical-splice
variants, against 3.87% for the liver enzymes ALT, AST and GGT. This pooled
contrast mixes trait class with GWAS ancestry and fine-mapping coverage: most
direct credible sets are European, most enzyme credible sets are non-European,
and the European enzyme GWAS were fine-mapped at 8 of 767 lead windows. Within
European studies the fractions were 29.88% and 6.92% (difference 0.2296, 95%
interval 0.0550 to 0.4236); within AMR, the one ancestry with comparable
fine-mapping coverage in both classes, they were 6.73% and 4.54% (difference
0.0219, -0.0489 to 0.1065). Neither stratum passed its calibration control, so
no trait-class call is made, and the pooled excess comes mainly from the three
PDFF GWAS. A trait-class comparison needs fine-mapping of all UKBB enzyme lead
windows and of the MVP European ALT, AST and NAFLD GWAS.""",True)
add('R14-03',AB,"""The strongest colocalization evidence for
423 of 436 genes came from a liver-enzyme trait and for only 13 from a direct
disease or liver-fat trait, so expression-based mapping of liver traits returns
predominantly a liver-enzyme regulatory map.""","""The strongest colocalization evidence for
423 of 436 genes came from a liver-enzyme trait and for 13 from a direct
disease or liver-fat trait; 22 genes have any direct-trait colocalization. This
split follows the number of signals tested: 171 direct and 3,844 enzyme signal
pairs were eligible, and 17.5% and 21.1% of them colocalized. The colocalized
map is dominated by liver-enzyme evidence because far more enzyme signals were
tested.""",True)
# intro
add('R02-01',IN,"""Fine-mapping the trait-associated regions shows that credible sets for direct
MASLD and liver-fat phenotypes carry 23.34% of their posterior probability in
protein-altering or canonical-splice variants, against 3.87% for the liver
enzymes ALT, AST and GGT (ratio 6.026671, 95% physical-cluster bootstrap interval
2.574593 to 13.510165). Separately, the strongest evidence for 423 of 436
colocalized genes comes from a liver-enzyme trait and for only 13 from a direct
disease or liver-fat trait, so the genetic map that expression-based
colocalization returns is largely a map of liver-enzyme regulation.""","""Pooled across ancestries, credible sets for direct MASLD and liver-fat
phenotypes carry 23.34% of their posterior probability in protein-altering or
canonical-splice variants, against 3.87% for the liver enzymes ALT, AST and GGT.
The pooled contrast mixes trait class with GWAS ancestry and fine-mapping
coverage; within European studies the fractions are 29.88% and 6.92%, within AMR
6.73% and 4.54%, and neither stratum supports a calibrated trait-class call.
Separately, the strongest evidence for 423 of 436 colocalized genes comes from a
liver-enzyme trait and for 13 from a direct disease or liver-fat trait (22 have
any direct-trait colocalization). Per eligible signal pair the colocalization
rate is similar (17.5% direct, 21.1% enzyme), so the enzyme majority reflects
the 22-fold larger number of enzyme signals tested.""",True)
# fig2 draft paragraph 101 (R02-03 scope sentence)
add('R02-03',F2,"We therefore quantify DNA architecture directly from all reliable Tier-1/2 SuSiE credible-set instances rather than",
"We therefore quantify DNA architecture directly from every reliable Tier-1/2 SuSiE credible set that fine-mapping produced rather than")
add('R02-03',F2,"The grouping limits repeated-study dominance and does not assert LD or statistical independence. Credible-set-weighted estimates remain a sensitivity analysis in the source table.",
"The grouping limits repeated-study dominance and does not assert LD or statistical independence. Credible-set-weighted estimates remain a sensitivity analysis in the source table. Fine-mapping coverage differs by stratum: every European direct-trait lead window has a credible set (47 of 47), but the European enzyme arm is 82 credible sets from UKBB GWAS Catalog lead windows (8 of 767 lead windows covered), and the MVP European ALT, AST and NAFLD GWAS were not fine-mapped.")
# fig2 paragraph 103 full replace
old103=open(D+F2).read().split('\n')[102]
assert old103.startswith('Taking protein-altering and canonical-splice')
new103=("Taking protein-altering and canonical-splice variants together, the pooled trait classes carry 23.3428% against 3.8732% of their posterior probability (difference 0.194695, 95% physical-cluster bootstrap interval 0.069229 to 0.339940; 233 physical 1-Mb clusters: 19 direct, 228 enzyme and 14 shared). This pooled value is not a trait-class contrast. Direct credible sets are mostly European (16 of 19 direct clusters contain a European set) and 6,126 of 6,208 enzyme sets are non-European, so trait class is confounded with GWAS ancestry, LD reference and fine-mapping coverage. Stratified by ancestry (Figure 2D, ancestry panel), European studies give 29.88% against 6.92% (difference 0.2296, 95% interval 0.0550 to 0.4236), AMR gives 6.73% against 4.54% (0.0219, -0.0489 to 0.1065), and all non-European studies pooled give 6.99% against 3.68% (0.0331, -0.0383 to 0.1274). The prespecified fixed-weight European-plus-AMR estimate is 0.1497 (0.0388 to 0.2718). None of these receives a trait-class call: the European enzyme arm was fine-mapped at only 8 of 767 lead windows, and the AMR and pooled estimates failed their random-split calibration controls. Within the direct class, the three European PDFF GWAS (37 credible sets, 12 clusters) carry 40.14% against 3.87% for all enzymes (difference 0.3626, 0.1835 to 0.5455), and NAFLD plus NASH (170 sets, 13 clusters) carry 12.44% (0.0857, -0.0164 to 0.2096). The direct excess is also concentrated: five of 19 direct clusters carry 76.0% of the summed direct cluster fractions, and the largest, one PDFF credible set at chr11:823,586, carries 22.5%; the PNPLA3 cluster, with 124 direct sets, has a mean of 0.117. An ungated sensitivity analysis matched on standardized lead-variant effect gives ratio 5.78 (95% interval 1.99 to 48.05; 27 matched enzyme loci), because direct loci come from smaller GWAS (median effective n 21,893 against 160,000) and carry larger standardized effects (median 0.0368 against 0.0183). Native credible-set member coordinates uniquely restore 715 sentinels whose decimal encoding had lost trailing zeros; historical September intervals from the faulty coordinates are superseded. The 57 adopted robustness comparisons (study, cluster, window and LD-block omissions) give ratios 4.368132 to 8.275691 but include no ancestry or trait stratification and are correlated sensitivities, not replications. No permutation p-value is reported. Comparing the trait classes within one ancestry requires SuSiE fine-mapping of all 767 UKBB ALT, AST and GGT lead windows and of the MVP European ALT, AST and NAFLD GWAS. COLOC r2 was separately adopted on October 14.")
add('R02-01',F2,old103,new103,True)
# fig2 line 105 (R14-03, R13-07)
add('R14-03',F2,"Separately, the strongest colocalization evidence for 423 of 436 genes comes from a liver-enzyme trait (ALT 137, AST 158, GGT 128) and for only 13 from a direct disease or liver-fat trait (NAFLD 8, PDFF 4, NASH 1). Because the table records one row per gene, this describes which trait class supplies each gene's best evidence rather than excluding weaker evidence from the other class. Expression-based colocalization applied to liver traits therefore returns predominantly liver-enzyme evidence.",
"Separately, the strongest colocalization evidence for 423 of 436 genes comes from a liver-enzyme trait (ALT 137, AST 158, GGT 128) and for 13 from a direct disease or liver-fat trait (NAFLD 8, PDFF 4, NASH 1). Because the table records one row per gene, this describes which trait class supplies each gene's best evidence; 22 genes have a direct-trait SuSiE PP.H4 above 0.5 (12 direct only and 10 with both classes; Figure 1C). The split follows the number of signals tested. Direct traits contributed 207 gated credible sets and 171 eligible signal pairs, enzyme traits 6,208 and 3,844; 30 of 171 (17.5%) and 810 of 3,844 (21.1%) pairs colocalized (Fisher odds ratio 0.80, nominal P = 0.29; pairs within a locus are not independent). Expression-based colocalization applied to liver traits therefore returns predominantly liver-enzyme evidence because many more enzyme signals were tested, not because direct-trait signals were observed to colocalize at a lower rate.")
# fig2 line 107 (R02-02, R01-10)
add('R02-02',F2,"The adopted coordinate-corrected cluster-weighted analysis also makes missing context visible: 36.49% of direct-trait and 77.85% of liver-enzyme noncoding posterior probability lacks resolved deposited regulatory context. Direct-trait ABC overlap is only 0.126%, and its bootstrap ratio interval reaches zero; sparse boundary behavior does not establish biological depletion. Consequently, the DNA result supports a broad noncoding architecture and an assay-coverage map, not a claim that every noncoding variant has a known regulatory mechanism.",
"Context annotation had failed for sites absent from the lookup table in non-European studies (0.7% of rows unresolved in EUR, 98.4% in AMR); re-annotating by GRCh38 position, with all 16,324 sites lifting 1:1, leaves no noncoding posterior probability without a context call. Weighting each credible set equally, as in the companion panel, 13.84% of direct-trait and 14.67% of liver-enzyme noncoding posterior probability lies within 2 kb of a TSS, 9.26% and 8.44% in liver-lineage open chromatin, 0.43% and 3.70% in a liver ABC enhancer, and 77.90% and 77.69% in sequence with none of these contexts. Within AMR, the one ancestry with comparable fine-mapping coverage in both classes, sequence with no deposited context showed no difference detected (direct minus enzyme -0.0957, 99% interval -0.344 to 0.126; Bonferroni over five contexts); promoter, ABC and accessible differences received no call because their split-sample calibration controls failed. Most noncoding fine-mapped probability therefore lies outside the deposited promoter, accessibility and ABC annotations, and the DNA result does not assign a regulatory mechanism to most noncoding signals.")
# ---------- legends
add('R02-03',LG,"""ancestries and six liver traits, followed by 360 trait-specific fine-mapped
loci (265 globally unique loci) and 996 COLOC genes. Brackets bind the""","""ancestries and six liver traits, with 360 trait-specific fine-mapped loci
(265 globally unique loci) and 996 COLOC genes counted separately. The COLOC
gene column is not derived from the fine-mapped loci: MVP European ALT, AST and
NAFLD were not fine-mapped, UKBB ALT, AST and GGT were fine-mapped at 6, 9 and
15 catalog lead variants, and 100 of the 436 multi-signal genes have a
colocalized lead within 250 kb of a fine-mapped window from the same study
(source table Fig2A_finemapped_support_by_cell). Brackets bind the""")
oldD=open(D+LG).read()
s=oldD.index('Coding plus splice mass is 23.34% versus 3.87%'); e=oldD.index('independent replications. No calibrated trait-class p-value is reported.')+len('independent replications. No calibrated trait-class p-value is reported.')
add('R02-01',LG,oldD[s:e],"""Pooled across ancestries, the proportion of fine-mapped posterior probability
carried by protein-altering or canonical-splice variants is 23.34% for direct
traits and 3.87% for enzyme traits (difference 0.1947, 95% physical-cluster
bootstrap interval 0.0692–0.3399; 10,000 draws, seed 20260817; 233 pooled
clusters: 19 direct, 228 enzyme, 14 shared). Because direct sets are mostly
European and enzyme sets mostly non-European, this pooled value is not a
trait-class contrast. **(D, ancestry panel)** Direct minus enzyme by stratum:
EUR 29.88% vs 6.92% (0.2296, 0.0550–0.4236), AMR 6.73% vs 4.54% (0.0219,
-0.0489–0.1065), fixed-weight EUR+AMR 0.1497 (0.0388–0.2718), all non-EUR
6.99% vs 3.68% (0.0331, -0.0383–0.1274). Thick bars are 95% and thin bars
98.333% intervals (Bonferroni over three prespecified members; 10,000 shared
draws, seed 20260923). No member receives a trait-class call: the EUR enzyme
arm was fine-mapped at 8 of 767 lead windows, and the AMR member failed its
random-split calibration control. **(D, direct-subclass panel)** Three EUR PDFF
GWAS (37 sets, 12 clusters) 40.14% and NAFLD plus NASH (170 sets, 13 clusters)
12.44%, each against all enzyme traits (3.87%); 95% intervals from 10,000
physical-cluster draws, seed 20260923. Five of 19 direct clusters carry 76.0% of
the summed direct cluster fractions. The 57-comparison robustness menu
(ratios 4.36813–8.27569) has no ancestry or trait stratification; these are
correlated sensitivities rather than independent replications. No calibrated
trait-class p-value is reported.""")
s=oldD.index('**(Noncoding-context companion, not assembled in canonical A-I)**'); e=oldD.index('gene, and regulatory-DNA context is distinct from transcript biotype.')+len('gene, and regulatory-DNA context is distinct from transcript biotype.')
add('R02-02',LG,oldD[s:e],"""**(Noncoding-context companion, not assembled in canonical A-I)** Deposited
annotations of noncoding posterior probability after re-annotating sites that
had failed lookup in non-European studies by GRCh38 position (16,324 sites, all
lifted 1:1); no noncoding posterior remains unresolved. Credible-set-weighted
fractions (direct versus enzyme): 13.84% versus 14.67% within 2 kb of a TSS,
9.26% versus 8.44% open in a liver lineage, 0.43% versus 3.70% overlapping a
liver ABC enhancer, and 77.90% versus 77.69% resolved with none of these
contexts. Context categories overlap and must not be summed. The earlier
displayed 52.11% and 60.43% unresolved fractions reflected failed lookup, not
missing annotation, and are withdrawn. Pooled values mix ancestries; within AMR
only the no-context difference passed calibration (no difference detected).
Regulatory-DNA context is distinct from transcript biotype.""")
add('R13-07',LG,"Coding plus splice mass is 23.34% versus 3.87%","(replaced by the R02-01 legend edit above; the phrase 'coding plus splice mass' is removed)")
# ---------- discussion
s0=open(D+DI).read(); s=s0.index('The fine-mapped architecture differs descriptively between trait classes.'); e=s0.index('(Rivas et al., 2015).\n\nExpression-based colocalization applied')+len('(Rivas et al., 2015).')
add('R02-01',DI,s0[s:e],"""Pooled across ancestries, protein-altering plus canonical-splice variants
carry 23.34% of fine-mapped posterior probability for direct MASLD/liver-fat
traits and 3.87% for enzyme traits. This pooled contrast is confounded with
GWAS ancestry and fine-mapping coverage. Within European studies the fractions
are 29.88% and 6.92%; within AMR they are 6.73% and 4.54%; neither stratum
supports a calibrated trait-class call, because the European enzyme GWAS were
fine-mapped at only 8 of 767 lead windows and the AMR estimate failed its
random-split calibration control. The pooled excess comes mainly from the three
PDFF GWAS and five of 19 direct clusters. Comparing the trait classes within
one ancestry needs fine-mapping of all UKBB enzyme lead windows and of the MVP
European ALT, AST and NAFLD GWAS. The earlier 6.69 ratio, its interval and
62-refit menu are historical and withdrawn after coordinate correction.
Separately, the strongest native colocalization evidence for 423 of 436 genes
comes from liver-enzyme traits and for 13 from direct disease/liver-fat traits;
per eligible signal pair the colocalization rate is similar (17.5% direct,
21.1% enzyme), so this split reflects how many signals of each class were
tested. Some protein-altering variants may act through protein function
without changing steady-state RNA, while others change RNA abundance through
nonsense-mediated decay or splicing (Rivas et al., 2015).""")
# ---------- METHODS
add('R02-01',ME,"""liver lineage-accessible chromatin were nonexclusive context annotations and
were not summed as probabilities.""","""liver lineage-accessible chromatin were nonexclusive context annotations and
were not summed as probabilities. Context was assigned by GRCh38 position after
1:1 liftOver of every member site, because lookup by study-specific variant key
failed for most non-European rows. Direct-minus-enzyme coding plus
canonical-splice differences were estimated pooled and within EUR and AMR
(fixed weights 16/26 and 10/26 for the combined estimate), with 10,000 shared
physical-cluster bootstrap draws (seed 20260923) and 98.333% intervals
(Bonferroni over three members). A stratum could carry a trait-class call only
if both classes had comparable fine-mapping coverage of their lead windows, the
minimum detectable difference was at most the pooled difference, and a
random-split control covered zero in at least 180 of 200 splits. Context
differences within AMR used 99% intervals (Bonferroni over five contexts).""")
# ---------- PAPER
add('R02-01',PA,"""- Direct MASLD and liver-fat credible sets have a larger descriptive
  protein-altering plus canonical-splice posterior fraction than liver-enzyme
  proxy credible sets: 23.3428% against 3.8732%, a 6.026671-fold ratio after""","""- Pooled across ancestries, direct MASLD and liver-fat credible sets have a
  larger descriptive protein-altering plus canonical-splice posterior fraction
  than liver-enzyme proxy credible sets. The pooled contrast is confounded with
  GWAS ancestry and fine-mapping coverage and is not a trait-class claim (H58,
  2026-10-15: EUR 29.88% vs 6.92%, AMR 6.73% vs 4.54%, no call in any
  stratum; pooled excess carried mainly by three PDFF GWAS). Pooled values:
  23.3428% against 3.8732%, a 6.026671-fold ratio after""",True)
add('R14-03',PA,"""  Separately, the colocalized map is predominantly a
  liver-enzyme map: the strongest evidence for 423 of 436 genes comes from a
  liver-enzyme trait and for 13 from a direct disease or liver-fat trait.""","""  Separately, the colocalized map is predominantly a
  liver-enzyme map because more enzyme signals were tested: the strongest
  evidence for 423 of 436 genes comes from a liver-enzyme trait and for 13 from
  a direct disease or liver-fat trait, 22 genes have any direct-trait
  colocalization, and per eligible signal pair the colocalization rate is 17.5%
  direct against 21.1% enzyme (171 and 3,844 pairs).""",True)
# ---------- STATUS
add('R02-01',ST,"""The descriptive direction persists, conditional on the studied GWAS portfolio.""","""The descriptive direction persists, conditional on the studied GWAS portfolio.
H58 (2026-10-15) found the pooled contrast confounded with ancestry and
fine-mapping coverage: EUR 29.88% vs 6.92%, AMR 6.73% vs 4.54%, no trait-class
call in any stratum; the 36.49%/77.85% (and 52.11%/60.43%) unresolved-context
values reflected failed lookup and are withdrawn (repaired unresolved 0%).""")
add('R02-02',ST,"""The companion retains instance-weighted unresolved fractions
52.11% and 60.43% and remains outside the assembled A–I panel set.""","""The installed companion still prints instance-weighted unresolved fractions
52.11% and 60.43%; H58 withdrew them (failed lookup, repaired value 0%), and a
corrected candidate companion is in
`figures/candidates/review-r3-fixes-20261015T112148Z/F2d/`. It remains outside
the assembled A–I panel set.""")
# ---------- RESULTS
add('R02-01',RE,"""fraction of 23.34% for direct traits and 3.87% for enzyme traits, ratio 6.03.""","""fraction of 23.34% for direct traits and 3.87% for enzyme traits, ratio 6.03,
pooled across ancestries. H58 (2026-10-15) shows this pooled value mixes trait
class with ancestry and fine-mapping coverage: EUR 29.88% vs 6.92% (difference
0.2296 [0.0550, 0.4236]), AMR 6.73% vs 4.54% (0.0219 [-0.0489, 0.1065]); no
stratum receives a trait-class call.""",True)
add('R02-02',RE,"""| Noncoding PIP mass with unresolved deposited context | 52.11% | 60.43% |""","""| Noncoding PIP mass with unresolved deposited context | 52.11% | 60.43% | (withdrawn 2026-10-15, H58: failed lookup; repaired 0% / 0%)""")
add('R02-02',RE,"""The five context differences (direct minus enzyme; 95% physical-cluster
intervals) are promoter 0.138659 [0.018337, 0.292354], ABC enhancer -0.005811
[-0.011888, -0.001089], accessible 0.059563 [-0.004648, 0.135635], resolved with
no context 0.226098 [0.076323, 0.379006], and unresolved context -0.413643
[-0.616094, -0.199384]. Direct ABC fraction is only 0.126% and its ratio
interval reaches zero; sparse/boundary behavior prevents a biological
depletion claim.""","""The five adopted context differences (promoter 0.138659, ABC -0.005811,
accessible 0.059563, resolved-no-context 0.226098, unresolved -0.413643) are
withdrawn: H58 (2026-10-15) showed the unresolved mass came entirely from failed
lookup in non-European studies. After position re-annotation (16,324 sites, 1:1
liftOver) the pooled differences are promoter 0.0995 [-0.028, 0.248], ABC
-0.0338 [-0.052, -0.017], accessible -0.0312 [-0.102, 0.046], resolved-no-context
-0.0896 [-0.229, 0.038], and unresolved 0. Pooled values mix ancestries; within
AMR (99% intervals, Bonferroni over five) only resolved-no-context passed
calibration (-0.0957 [-0.344, 0.126], no difference detected).""")
add('R02-02',RE,"""remains outside A-I and displays credible-set-weighted fractions: unresolved
52.11% / 60.43%, promoter 7.87% / 4.46%, accessible 7.32% / 4.15%, ABC
0.23% / 1.57% (direct / enzyme).""","""remains outside A-I and displays credible-set-weighted fractions: unresolved
52.11% / 60.43%, promoter 7.87% / 4.46%, accessible 7.32% / 4.15%, ABC
0.23% / 1.57% (direct / enzyme). H58 repaired values (candidate panel): unresolved
0% / 0%, promoter 13.84% / 14.67%, accessible 9.26% / 8.44%, ABC 0.43% / 3.70%,
resolved-no-context 77.90% / 77.69%.""")
add('R14-06',RE,"""| Smallest difference-interval lower bound | **0.022373** | All 57 difference intervals above zero; each has 10,000 finite draws; these are correlated sensitivities, not 57 independent replications |""","""| Smallest difference-interval lower bound | **0.022373** | All 57 difference intervals above zero; each has 10,000 finite draws; these are correlated sensitivities, not 57 independent replications |
| Direct-cluster concentration (2026-10-15) | top 5 of 19 clusters carry 76.0% of summed direct cluster fractions; chr11:823,586 (one PDFF set) 22.5% | Ratio after removing top 1-5 clusters: 4.93, 4.15, 3.45, 2.69, 1.96 (biased diagnostic, not a test) |
| Direct subclass (2026-10-15) | PDFF 40.14% (diff 0.3626 [0.1835, 0.5455]); NAFLD+NASH 12.44% (0.0857 [-0.0164, 0.2096]) vs enzyme 3.87% | 10,000 physical-cluster draws, seed 20260923 |""")
# verify
bad=[]
for x in E:
    if x['fid']=='R13-07' and x['old_text'].startswith('Coding plus'): 
        x['note']='subsumed'; continue
    t=open(D+x['file']).read(); n=t.count(x['old_text'])
    if n!=1: bad.append((x['fid'],x['file'],n,x['old_text'][:60]))
    else: x['line']=str(t[:t.index(x['old_text'])].count('\n')+1)
print('bad',bad); print(len(E))
json.dump(E,open('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F2d/text_edits.json','w'),indent=1)
