import re
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
E=[
('R06-06','docs/manuscript/05_figure_legends.md',384,"with mRNA; 2,843 are not protein-significant. Spearman rho is 0.457835 and"),
('R06-06','docs/manuscript/05_figure_legends.md',385,"direction agreement is 65.9136%, rounded to 66% on the panel."),
('R06-06','docs/manuscript/draft/fig5.md',13,"liver proteome (58 patients). Of the genes reaching protein significance, most\nchange in the same direction as their transcript; nonsignificance is not\ninterpreted as falsification."),
('R06-06','docs/RESULTS.md',3064,"65.9136%, rounded to 66% on the panel. OLS slope=0.455025,"),
('R06-08','docs/manuscript/05_figure_legends.md',403,"MASLD/liver-fat pairs (n=28) from liver-enzyme pairs (n=739) in the adopted\nOctober 14 r2 PDF. The historical panel had 29 and 787 pairs."),
('R06-08','docs/manuscript/draft/fig5.md',35,"Figure 5D intersects the regulatory posterior of 767 colocalized gene-study\npairs with open chromatin from two liver snATAC sources across five shared\nlineages, keeping direct MASLD/liver-fat and liver-enzyme traits separate. Peak\noverlap localizes regulatory DNA; it does not assign a variant to a neighboring\ngene, and absent overlap is unresolved context rather than a negative test."),
('R06-08','docs/RESULTS.md',2940,"28 direct MASLD/liver-fat pairs from 739 liver-enzyme pairs (767 total), read"),
('R06-08','docs/STATUS.md',829,"5D uses 767 r2 gene-study pairs (28 direct, 739 enzyme)."),
('R06-08','figures/main/fig5_molecular_context/README.md',28,"colocalized gene-study pairs (436 genes; 28 direct MASLD/liver-fat and 739\nliver-enzyme pairs), then sums normalized SNP.PP.H4"),
('R06-08','figures/main/fig5_molecular_context/README.md',33,"traits remain separate, and\nboth-source, one-source, partial, no-overlap, and untestable states remain\nexplicit. Peak overlap localizes regulatory DNA; it does not assign a target\ngene or establish functional regulation."),
('R11-02','docs/manuscript/05_figure_legends.md',405,"both-source, one-source, partial, and no-overlap states; all pairs met the 0.95"),
('R11-02','docs/manuscript/05_figure_legends.md',407,"source table. Peak overlap localizes regulatory DNA but does not assign a\ntarget gene or establish functional regulation."),
('R11-02','docs/RESULTS.md',2946,"are partial rather than both-source accessible. Peak overlap is regulatory-DNA"),
('R11-02','docs/manuscript/METHODS.md',1121,"separate. Partial or absent peak overlap is unresolved context, not a negative\nfunctional test."),
('R06-09','docs/manuscript/draft/fig5_discussion.md',36,"GNMT, MAT1A, and CYP2C19\ncarry concordant disease-state RNA, cross-sectional stage, and tissue-protein\nevidence, routing the hepatocyte metabolic-function hypothesis to\ncontext-matched perturbation and rescue."),
('R11-07','docs/manuscript/05_figure_legends.md',379,"caption. GNMT, MAT1A, and CYP2C19 are lower at both RNA and protein levels and\nare interpreted as severity-linked reduction in hepatocyte metabolic gene and\nprotein abundance."),
('R11-07','docs/manuscript/draft/fig5.md',23,"as pooled RNA. This is assay-native support for a severity-linked reduction\nin hepatocyte metabolic gene and protein abundance, not evidence that the three genes share one\ninherited regulatory mechanism."),
('R11-07','docs/RESULTS.md',553,"### Loss of hepatocyte metabolic function: GNMT–MAT1A–CYP2C19"),
('R11-07','docs/RESULTS.md',576,"a recurring example of reduced abundance of hepatocyte metabolic genes and proteins in\nestablished disease. It does not establish"),
('R11-07','docs/manuscript/05_figure_legends.md',399,"rows are descriptive rather than independent validation."),
]
for fid,f,ln,old in E:
    t=open(B+f).read(); n=t.count(old)
    if n:
        i=t.index(old); L=t[:i].count('\n')+1
    else: L=None
    print(fid,f,ln,'count',n,'found_at_line',L)
