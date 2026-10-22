EDITS = [
("R04-01","docs/manuscript/draft/abstract.md","48-51","""Donor-aware single-cell analysis and a prespecified 117-program registry
resolved remodeling into multicellular programs while ambient-RNA sensitivity
limited lineage assignment. After source-diagnosis correction, Ductular injury
alone met the primary 117-program rule; no program passed HC3 family correction."""),
("R04-01","docs/manuscript/draft/fig4.md","68-78","""Figure 4C shows the three positive and three negative program effects with the
smallest complete-family BH q per lineage, with all 117 programs retained in
the supplementary view. ECM/IGFBP7 and ductular-injury/BICC1, both discovered in
hepatocyte profiles, are among the 37 whose content contradicts that label.
After correcting two GSE174748 explants to source-documented Cirrhosis and
applying the unchanged cirrhosis exclusion, ordinal-stage models use at most
62 biological donors across five datasets with dataset as a covariate.
Only BICC1 passes the primary 117-program BH/stability rule (beta=0.4021,
q=0.03763); ECM does not (beta=0.3515, q=0.08546). Neither passes HC3 family
correction (BICC1 q=0.1558; ECM q=0.1652). The previous 64-donor two-program
selection is historical."""),
("R04-01","docs/PAPER.md","266-269","""On at most 62 source-reported donors,
   only Ductular injury passes the primary 117-program BH/stability rule
   (β=0.402086, q=0.037633); Stromal ECM does not (β=0.351525, q=0.085458).
   Neither passes HC3 family correction."""),
("R04-01","docs/RESULTS.md","627-634","""The source-diagnosis correction completed October 14 uses at most 62
source-reported donors across five datasets, tests all 117 programs in one BH
family, and adjusts for dataset. The two GSE174748 disease explants are
cirrhotic in the original [Filliol study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9949942/),
so the unchanged three-stage rule excludes them. Stromal ECM has disease
β=0.351525, primary q=0.085458 and HC3 q=0.165240; it no longer passes primary
selection. Ductular injury has β=0.402086, primary q=0.037633 and HC3 q=0.155826;
it is the sole primary-selected program. Neither passes HC3 family correction."""),
("R04-01","docs/RESULTS.md","1600-1601","""| Programs meeting primary association/stability criteria | 1 | Diagnosis-qualified 62-donor maximum on the ambient-uncorrected atlas; BICC1 only; BH over all 117 programs |
| Programs passing HC3 family correction | 0 | Corrected ECM HC3 q=0.165240 and BICC1 q=0.155826 |"""),
("R04-01","docs/RESULTS.md","1638-1644","""| Atlas label | Program | Primary q | HC3 q | Source-diagnosis-refit state |
|---|---|---:|---:|---|
| Hepatocyte | Stromal ECM | 0.085458 | 0.165240 | Not primary-selected |
| Hepatocyte | Ductular injury | 0.037633 | 0.155826 | Primary-selected; fails HC3 family correction |

These two rows are the October 14 diagnosis-qualified, 62-donor maximum
single-cell refit."""),
("R04-01","docs/RESULTS.md","1767-1768","""lncRNA-independent. Both source-diagnosis-refit primary-selected programs passed this
lncRNA-removal sensitivity:"""),
("R04-01","docs/RESULTS.md","2785-2787","""117-program stage context has one primary-selected program, BICC1/Ductular
injury, zero HC3-selected programs, 116 ordinary nonsignificant programs and
zero strict negative calls."""),
("R04-01","docs/STATUS.md","22-25","""disease explants to Cirrhosis. The unchanged three-stage rule excludes them:
maximum primary n=62, BICC1 alone primary-selected, zero HC3-selected;
maximum ambient common n=60. Figure 4C and the full 117-program supplement
use the corrected effects."""),
("R04-01","docs/STATUS.md","822","""primary maximum n=62, BICC1 alone selected (q=0.037633), ECM not selected (q=0.085458), zero HC3-selected."""),
("R04-01","docs/README.md","12-13","""Original-source diagnosis correction also refits the fixed 117-program
stage models: maximum n=62, BICC1 alone primary-selected, zero HC3-selected."""),
("R04-01","docs/ROADMAP.md","19-21","""The fixed-score refit has maximum primary n=62,
BICC1 alone selected and zero HC3-selected; its ambient common universe has
maximum n=60."""),
("R04-01","docs/manuscript/README.md","28","""Diagnosis-qualified 62/60-donor primary/ambient fits; Ductular injury alone primary-selected, zero HC3-selected; fixed weights and unresolved lineage origin"""),
("R04-01","docs/manuscript/draft/fig5_discussion.md","44-46","""The diagnosis-qualified refit retains primary support only for Ductular injury;
Stromal ECM no longer passes the 117-program BH/stability rule, and neither
passes HC3 family correction."""),
("R04-01","docs/manuscript/draft/01_intro.md","99-101","""Source-diagnosis correction retains only Ductular injury
under the primary 117-program rule; neither focal program passes HC3 family
correction."""),
("R04-02","docs/manuscript/METHODS.md","651-652","""The corrected primary ordinal models therefore use at most 62 donors across
five datasets, adjusting for dataset; the former 64-donor fit is historical."""),
("R04-03","docs/RESULTS.md","1606-1608","""The corrected disease-state refit contains at most 62
biological donors across GSE185477, GSE202379, GSE244832, GSE174748, and the
two analyzed GSE189600 donors, with no Liver Atlas rows."""),
("R04-04","docs/manuscript/draft/fig4.md","79-84","""ECM raw beta=0.3509, corrected beta=0.1279, and
paired delta=-0.2230 (delta BH q=0.4781); BICC1 raw beta=0.3846,
corrected beta=0.4010, and delta=0.01642 (delta BH q=1). The former
62-donor ambient values are historical. BICC1 remains
`attenuation_indeterminate`; ECM is
`not_applicable_not_frozen_disease_association` after losing primary selection."""),
("R04-04","docs/PAPER.md","271-273","""Neither paired change passes the 117-program
   BH family; ECM's selection-conditioned attenuation classification is no
   longer applicable."""),
("R04-04","docs/RESULTS.md","665-666","""Its paired correction-induced change is β=-0.223006
(BH q=0.478071; HC3 q=0.337646), and its"""),
("R04-04","docs/RESULTS.md","671-673","""BICC1 remains
`attenuation_indeterminate`; ECM's selection-conditioned attenuation state is
not applicable after losing primary selection."""),
("R04-05","docs/manuscript/05_figure_legends.md","228-238","""**(C)** The three positive and three
negative Hotspot effects with the smallest complete-family BH q in each of the
four non-T lineages, shown across Healthy, Steatosis, and Steatohepatitis and
paired with signed disease-stage β estimates and 95% confidence intervals. T
cell programs are omitted from this compact main panel; the complete 117-program
registry view, including T cell programs, is shown in Figure S4. The corrected
62-donor primary fit excludes the two GSE174748 cirrhotic explants under the
original stage rule. Only Ductular injury passes the 117-program BH/stability
rule (beta 0.4021, q=0.03763); Stromal ECM does not (beta 0.3515, q=0.08546).
Neither passes HC3 family correction (q=0.1558 and 0.1652, respectively).
The prior 64-donor two-program selection is historical."""),
("R04-06","docs/manuscript/05_figure_legends.md","296-298","""The native forest
shows cohort and cohort-omission disease-grade effects in hepatocytes; its
intervals are nominal sensitivity intervals, not new complete-family calls."""),
("R04-07","docs/manuscript/draft/fig4.md","85-86","""Taken with the content audit, this prevents any
hepatocyte-intrinsic claim from the original atlas."""),
("R04-07","docs/RESULTS.md","686-687","""These results support same-atlas descriptive localization, not transcript
origin or cell autonomy."""),
("R04-07","docs/manuscript/METHODS.md","686-687","""reported. This analysis describes cellular context within the same atlas and
does not test transcript origin, secretion, or cell autonomy."""),
("R04-08","docs/manuscript/draft/fig4.md","107-111","""Figure 4D asks the separate question of cellular context. In donor-paired
same-atlas contrasts, Stromal ECM/IGFBP7 is highest in fibroblasts and Ductular
injury/BICC1 is highest in cholangiocytes. The ten prespecified contrasts pass
both ordinary and HC3 family criteria for those focal lineages. This is
descriptive localization, not transcript origin or cell autonomy."""),
("R04-09","docs/manuscript/draft/fig4.md","54-56","""Scored against donor-level per-lineage expression, **80 of 117
programs match their declared lineage and 37 do not**, and the mismatch is
systematic:"""),
("R04-09","docs/manuscript/05_figure_legends.md","245-246","""Eighty match their discovery lineage and
37 do not."""),
("R04-10","docs/manuscript/draft/01_intro.md","97-99","""ECM/IGFBP7 and ductular-injury/BICC1 test whether a disease-associated program,
its cell-of-origin assignment, and its intact-tissue organization support the
same conclusion."""),
("R10-01","docs/manuscript/draft/fig4.md","193-195","""Likewise, the existing Hotspot matched-background comparison
does not establish exceptional disease specificity (fibrosis 26/114 versus
null mean 18.445, nominal P=0.134; NAS 4 versus 4.225, P=0.512)."""),
("R10-01","docs/RESULTS.md","2146-2147","""the vocabulary null (26 of 114
fibrosis-supported against 18.4 expected, p = 0.134) now has a cause,"""),
("R10-01","docs/RESULTS.md","1635-1636","""The 27/113 belongs to the 469-record working roster;
the 774-unit arms are in the Figure 4E paragraph below."""),
("R13-04","docs/manuscript/draft/01_intro.md","54-56","""Donor-aware
single-cell analyses and prespecified gene programs resolve multicellular
remodeling, while protein,"""),
("R13-04","docs/manuscript/draft/fig1.md","44-45","""Figure 4
resolves donor-level multicellular programs and regulatory activity."""),
("R13-04","docs/PAPER.md","122-123","""2. Multicohort bulk and donor-aware single-cell analyses resolve established
   disease as multicellular, cross-sectional remodeling."""),
("R13-04","docs/manuscript/draft/title.md","10-11","""**A multimodal human MASLD Resource resolves the architecture of inherited
risk and the coordination of multicellular remodeling**"""),
("G1-01","docs/manuscript/METHODS.md","27-28","""The 50-cell gate is primary; 20 cells and
annotation confidence ≥0.90 are sensitivities."""),
("G1-02","docs/manuscript/METHODS.md","636-637","""after two prespecified low-cell
GSE189600 libraries are excluded."""),
("G1-03","docs/RESULTS.md","1599","""| Cell types represented | 5 | Hepatocyte, macrophage, fibroblast, cholangiocyte, T cell |"""),
("G1-04","docs/manuscript/METHODS.md","637-638","""Dataset identity was retained during
integration."""),
]
