<!-- Read-only assembled copy of release 2026-10-20-r8, built 2026-10-20 from docs/manuscript/draft/*.md, 05_figure_legends.md and METHODS.md. Edit the section files, not this file. -->

# Title page

**A multimodal human MASLD Resource records liver-trait genetic evidence and
disease-state liver expression side by side for each gene**

Authors: James Lee^1,2,\*^, Victor Passanisi^1,2,\*^, Neha Saravanan^1,2^, Neville E. Sanjana^1,2,\#^

Affiliations: ^1^New York Genome Center, New York, NY, USA. ^2^Department of Biology, New York University, New York, NY, USA.

\* These authors contributed equally.

Lead contact and corresponding author: Neville E. Sanjana (\#), nsanjana@nygenome.org

Keywords: [to be supplied]

# Summary

Gene-level evidence from inherited variation and from the disease-state liver transcriptome of established metabolic dysfunction-associated steatotic liver disease (MASLD) is rarely recorded side by side with what each assay could test. Here, we assembled a human MASLD Resource combining bulk liver RNA sequencing from five cohorts, analyzed as 774 genetic-identity clusters rather than verified individuals, with liver eQTL colocalization of 35 liver-trait GWAS strata. We identified 1,854 disease-state calls and 436 multi-signal colocalized genes. For 423 of the 436, the strongest colocalization came from a liver-enzyme trait rather than a direct disease or liver-fat trait. Within expression deciles, we detected neither enrichment nor depletion of genes in both sets, although expression was not fully controlled. Most disease-state calls without colocalization were untested, not genetic negatives. The MASLD Gene Catalog records genetic and disease-state evidence separately for 28,233 gene identities, distinguishes untested from unsupported, and assigns each gene a next experiment.

# Introduction

Metabolic dysfunction-associated steatotic liver disease (MASLD) is a common chronic liver disorder that spans steatosis, metabolic dysfunction-associated steatohepatitis (MASH), fibrosis and cirrhosis [ref: MASLD and MASH nomenclature consensus]. Fibrosis stage is strongly associated with liver-related and overall outcomes [ref: fibrosis stage and clinical outcomes in MASLD cohorts]. Yet how inherited genetic variation relates to the multicellular remodeling of established human disease is still incompletely mapped. Liver samples with similar histological diagnoses can contain different mixtures of hepatocyte stress, immune activation, stromal remodeling and vascular dysfunction, so tissue and cell context are needed to interpret them. Histology is therefore an essential clinical anchor, but it does not fully describe molecular state. A data-driven molecular continuum, following the framework of Kamzolas et al. (2026; doi:10.1038/s42255-026-01543-7), can add resolution within recorded stages without replacing pathological assessment.

Human genetics and disease-tissue transcriptomics approach this biology from different directions. Genome-wide association studies (GWAS) and fine-mapping identify inherited trait-associated variation. This variation can include protein-altering variants, noncoding variants with regulatory evidence and signals whose molecular mechanism is unresolved [ref: statistical fine-mapping methods]. Colocalization with molecular quantitative trait loci (QTL) can connect some trait associations to gene regulation [ref: colocalization methods, including multi-signal colocalization]. However, colocalization alone does not establish which gene mediates a trait association. It also does not require the gene to stay differentially expressed once disease is established. Some protein-altering variants may act through protein function without changing steady-state RNA, while others change RNA abundance through nonsense-mediated decay or splicing (Rivas et al., 2015). Conversely, differential expression captures active changes in tissue state but cannot tell inherited drivers apart from adaptive responses, shifts in cell composition or consequences of chronic injury. The two maps therefore answer different but connected biological questions.

Recent studies have produced large liver expression QTL (eQTL) resources [ref: human liver eQTL resources], multicohort transcriptomic analyses [ref: multicohort MASLD liver transcriptome analyses], single-cell and spatial MASLD atlases [ref: single-cell and spatial human MASLD liver atlases] and multilayer gene-prioritization frameworks [ref: multilayer MASLD gene-prioritization frameworks]. These resources have nominated risk loci, disease-associated cell states, fibrosis programs and potential therapeutic mechanisms. Kendall et al. (2023) linked genes to clinical outcomes in a multimodal MASLD database, and LiRNA (Alvarez-Sola et al., 2026) made harmonized human liver RNA-seq interactively queryable. Neither carries a fine-mapped inherited-variation layer next to its expression layer, so neither can say whether the two point to the same genes. Other studies have linked parts of the two maps. A random-effects meta-analysis of 10 biopsy datasets (1,058 samples) integrated MASLD and MASH differential expression with GWAS and coexpression modules (Piras and DiStefano, 2024). A GWAS meta-analysis paired with liver RNA-seq and genotypes from 504 participants nominated the long noncoding RNA (lncRNA) TRIB1AL by colocalization (Gobeil et al., 2026a). Surrogate-trait GWAS, colocalization and CRISPR interference (CRISPRi) in hepatocytes identified VKORC1, TNKS, LYPLAL1 and GPAM as regulators of lipid accumulation (Saliba-Gustafsson et al., 2025). Partitioned polygenic scores separated liver-confined from systemic MASLD (Jamialahmadi et al., 2024).

None of these studies keeps fine-mapped multi-ancestry colocalization, multicohort disease-state expression and donor-aware cell and tissue context side by side for each gene. Here, donor-aware means that the donor, not the cell, is the unit of inference. What is missing is a per-gene record of both kinds of evidence that also states, for each layer, whether the gene was untested or was tested without support. Such a record has to keep the biological meaning, coverage, source and unit of inference of each assay, rather than reduce them to a count of modalities.

Here, we assembled a multimodal human MASLD Resource that keeps these two kinds of evidence side by side (Figure 1). The transcriptomic arm includes bulk liver RNA sequencing from nine cohorts. Of these, five cohorts have controls, and their 774 genetic-identity clusters entered a pooled, cohort-adjusted disease-versus-control analysis (Figure 3). These clusters are not verified independent people. The genetic arm combines prespecified GWAS strata across ancestry groups with liver eQTL colocalization, fine-mapped coding and noncoding context, and a reanalysis of a published liver chromatin-accessibility QTL map (Currin et al., 2025). It keeps direct disease and imaging phenotypes separate from liver-enzyme proxy traits (Figure 2). Single-cell analyses score prespecified gene programs in each cell lineage, with donors as the unit of inference (Figure 4). Protein, chromatin-accessibility and spatial assays add molecular and physical context, each in its native unit (Figure 5). We retain fibrosis stage in the primary severity analyses. We also add an externally anchored transcriptomic continuum to describe molecular variation among participants assigned the same stage.

We follow two recurring examples through these layers. Both were selected from prespecified analyses, not from a count of supporting modalities. GNMT has a single-signal colocalization with NAFLD, a direct disease trait, but no multi-signal colocalization call, and its NAFLD association is sub-threshold (p = 2.6 × 10⁻⁴). Its Catalog genetics call is therefore indeterminate, and its local regulatory haplotype is unresolved. Together with MAT1A and CYP2C19, GNMT has lower disease-state RNA, lower abundance at higher fibrosis stage in cross-sectional comparisons, and lower liver protein. We do not claim that this three-gene metabolic state reflects one inherited mechanism. The second example is a pair of programs, ECM/IGFBP7 and ductular-injury/BICC1. We use them to test whether a program's disease association, its cell-of-origin assignment and its organization in intact tissue support the same conclusion. In bulk RNA from the 364 units, controls included, of the three cohorts that record both an exact fibrosis stage and NAS, the fibrosis association of ECM/IGFBP7 does not reach Benjamini-Hochberg (BH) q < 0.05. In single-cell data from 60 donors, neither program has a cross-sectional stage association that passes BH correction across all 117 programs, with p values from stage labels shuffled within each dataset (Supplementary Note 2). Ambient RNA is RNA from other cells that contaminates single-cell profiles. Program scores before and after ambient correction do not let us assign either program to hepatocytes. In spatial transcriptomics, only ECM/IGFBP7 is more spatially clustered than matched control genes. Ductular-injury/BICC1 is indeterminate (Figure 5F).

We then compared the two maps gene by gene. Of the 14,704 genes that both analyses list, multi-signal colocalization evaluated 3,769. Among these 3,769, 426 were colocalized and 220 were disease-state calls. The two sets shared 26 genes, against 21.2 expected within expression deciles (Mantel-Haenszel odds ratio 1.29, 95% CI 0.84 to 2.01). Because expression was not fully controlled, we state no direction. The overlap is also not informative on its own. The interval includes 1, the unstratified test was underpowered across a wide range of odds ratios, and two sets of these sizes are mostly disjoint whether or not the underlying biology differs. Across all 14,704 shared genes, 935 of the 1,129 disease-state calls without colocalization had a single-signal colocalization result but no multi-signal result, so they are not genetic negatives (Figure 1; Supplementary Note 3). The Gene Catalog gives such genes an indeterminate genetic call. Its untested class holds genes untestable in either domain, including 698 of the 699 disease-state calls that lie outside the 14,704 genes and have no tier-1/2 colocalization record (Figure 6D). Matched for expression, colocalized genes are not more hepatocyte-specific than disease-state genes, so the prespecified directional prediction failed. The difference we observed ran the other way, but that comparison is post hoc and does not establish confirmatory class separation (Supplementary Note 3). We detected no difference between the two classes in the liver protein response to fibrosis stage (73 pairs). This is not evidence of equivalence. Pooled across ancestries, fine-mapped credible sets for direct MASLD and liver-fat phenotypes carry 23.34% of their posterior probability in protein-altering or canonical-splice variants. The value is 3.87% for the liver enzymes alanine aminotransferase (ALT), aspartate aminotransferase (AST) and gamma-glutamyltransferase (GGT). This contrast mixes trait class with GWAS ancestry and fine-mapping coverage, so it receives no trait-class call (Figure 2; Supplementary Note 4).

The resulting MASLD Gene Catalog gives each protein-coding gene, lncRNA transcript and unresolved locus separate evidence states for inherited variation and for disease-state RNA. It links genes to cell programs by membership, attaches protein, chromatin and spatial context at the program level, and assigns the next discriminating experiment (Figure 6). The record separates an inherited candidate for perturbation from a disease-state change that could serve as its readout, and names the assay needed to connect the two. Finally, we score predictive models against each task's own baseline and report each score beside a per-study record of the model's training exposure (Figure 7). This per-gene record, which keeps the two kinds of evidence apart and separates untested from unsupported, is the standalone contribution of the Resource.

# Results

## A multimodal human MASLD Resource records inherited variation and disease remodeling side by side

To compare inherited genetic evidence with disease-state evidence gene by gene, we assembled a multimodal human MASLD Resource. It records trait-associated inherited variation, multicohort disease remodeling, and donor-aware cell and tissue context side by side (Figure 1A). For each result, the Resource lists the source cohort, assay, phenotype, biological unit and set of genes tested, whether the same data entered another analysis or source study, the evidence state, the main remaining alternative explanation, and the next discriminating experiment (Figure 1A). The transcriptomic data comprise 1,257 libraries that passed technical quality control, from nine cohorts (GSE126848, GSE130970, GSE135251, GSE162694, GSE167523, GSE174478, GSE193066, GSE213621 and GSE240729). The pooled disease-versus-control model uses 774 genetic-identity clusters from the five cohorts that include controls (138 control and 636 disease clusters). These clusters are not verified independent donors. The genetic data span prespecified liver GWAS strata from several ancestry groups. We kept direct MASLD and liver-imaging traits separate from the liver enzymes ALT, AST, and GGT, which we treat as proxy phenotypes. Donor-collapsed single-cell data, prespecified gene programs, and liver protein, chromatin-accessibility, and spatial measurements add context in each assay's native units. They are not extra votes in a single target score. Throughout, we define replication by participants or biological donors; cells, spots, sequencing runs, guides, and arrays do not add to the biological sample size.

We first asked how often the two maps name the same gene. We intersected the two gene sets within one shared testing universe: the 3,769 genes (stable Ensembl identifiers) evaluated by multi-signal colocalization (COLOC run on SuSiE fine-mapped signals, which allows more than one association signal per locus; Figure 1B). Of these, 426 were colocalized with a liver eQTL (PP.H4 > 0.5, where PP.H4 is the posterior probability that the GWAS and eQTL associations share one variant), 220 were canonical disease-state genes (BH q < 0.05 and |log2FC| > 0.5 in the pooled disease-versus-control fit), and 26 were in both. Within expression deciles, 21.2 shared genes were expected (Mantel-Haenszel OR 1.295, 95% CI 0.835 to 2.007; exact stratified *p* = 0.2755, one prespecified test). The interval includes 1. The unstratified test had less than 80% power to detect odds ratios between about 0.46 and 1.77. As a control, a prespecified list of non-colocalized genes matched to the colocalized genes on expression carried 13 disease-state genes against 21.3 expected within the same deciles. Expression level is therefore not fully controlled in this comparison, and we read no direction from it.

Several groups of genes fall outside this comparison. Another 935 canonical disease-state genes had no multi-signal COLOC result; we do not count them as genetic negatives. Seven genes are SuSiE-untestable because, in one tier-1/2 study (tier 1, direct NAFLD, NASH and liver-fat GWAS; tier 2, liver-enzyme GWAS), no GWAS-eQTL signal pair passed the shared-posterior check. This check requires both signals in a pair to keep at least half of their posterior probability on variants that the GWAS and eQTL fits share. These seven genes are not counted as positive and are also outside the 3,769 genes. Their exclusion does not establish a negative SuSiE result. The 3,769 genes are drawn from the 14,704 genes, matched by stable Ensembl identifier, that have both a disease-versus-control bulk test and a liver-eQTL record (Methods). The Gene Catalog counts genes differently. Grouping by gene symbol gives 14,720 genes, and setting aside 31 identities with unresolved annotation and seven genetically untestable records leaves 14,682 Catalog genes. These gene totals are not interchangeable. Supplementary Note 5 reports a signed comparison, made in a separate universe of 14,703 genes, between genetically predicted liver expression for ALT and the disease-versus-control statistic [ref: OTTERS TWAS method].

We next asked which liver traits supply the colocalization support (Figure 1C). We counted support from any eligible SuSiE gene-trait pair, not only from the strongest trait. Of the 436 primary SuSiE genes (426 of which fall in the 3,769-gene comparison), 414 were supported only by liver-enzyme GWAS, 12 only by direct-disease or liver-fat GWAS, and 10 by both. The two trait classes were not tested equally. Of the 414 enzyme-only genes, 46 were SuSiE-tested in a NAFLD, NASH, or MRI proton density fat fraction (PDFF) GWAS and 368 in none. Of the 12 direct-only genes, 11 were SuSiE-tested in an enzyme GWAS. Among the 362 genes SuSiE-tested in both trait classes, 56 were enzyme-supported and 21 direct-supported, including 10 supported by both. We describe the fine-mapped architecture of each trait class in Figure 2.

The Resource is organized as a graph of source dependencies and an analysis workflow that ends in a decision for each gene (Figure 1D). Figure 2 maps phenotype- and ancestry-aware inherited variation and its candidate molecular links. Figure 3 defines multicohort, cross-sectional disease remodeling and where it meets genetics. Figure 4 reports donor-level program scores in five cell lineages with their bulk-tissue context; transcription-factor activity is in the supplement. Figure 5 places candidate genes and programs in molecular and physical tissue context. Figure 6 organizes the two maps in the MASLD Gene Catalog and names the smallest experiment that could resolve each remaining alternative. Figure 7 asks what predictive models learn from the Resource. No step requires agreement across all modalities, and no aggregate rank is used as a calibrated measure of biological truth.

Two sets of examples recur through these figures. GNMT, MAT1A, and CYP2C19 connect the RNA, cross-sectional stage, and protein views of hepatocyte metabolic remodeling. The ECM/IGFBP7 and ductular-injury/BICC1 programs connect prespecified donor-level programs to spatial interpretation against matched nulls. We report the colocalization state of each example, but that state does not determine its role in disease remodeling.

Taken together, this design lets a reader trace why a gene was nominated, see what has been tested for it, and choose an experiment that could distinguish the remaining alternatives. None of these steps depends on the outcome of a perturbation screen. A separate prospective perturbation study may use the Resource, but its gene list, guides, phenotypes, hits, and mechanisms are not part of this paper.

## Multi-signal liver eQTL colocalization links 436 genes to 35 liver-trait GWAS strata

To map inherited variation relevant to MASLD, we assembled a prespecified set of direct-disease, MRI-derived liver-fat, and liver-enzyme GWAS from several ancestry groups (Figure 2A) [ref: source publications for the 35 GWAS strata]. We kept direct disease and imaging traits apart from the ALT, AST, and GGT proxy phenotypes. Enzyme associations report liver biology but are not equivalent to direct MASLD susceptibility. The registry contains 35 tier-1/2 strata from five GWAS ancestries and six liver traits. Fine-mapping retained 360 trait-specific loci, which correspond to 265 globally unique loci [ref: SuSiE fine-mapping].

Colocalization against the Broadaway liver-eQTL panel (Broadaway et al., 2024) gave two maps that differ in what they can observe (Figure 2B). As in Figure 1, each result records phenotype applicability, eQTL and LD coverage, method eligibility, evidence state, unresolved alternative, and next discriminating experiment. Both methods report PP.H4, the posterior probability that the GWAS and eQTL associations share one variant. Single-variant ABF and multi-signal SuSiE colocalization have different eligibility requirements [ref: coloc ABF method; coloc-SuSiE method]. SuSiE needs converged eQTL and GWAS fits and usable LD, so a missing multi-signal result marks limited observability, not a negative result.

The genetic union contains 996 genes: 436 primary multi-signal genes and 560 single-signal-only genes. Seven genes remain SuSiE-untestable, three of them in the single-signal-only set, and they are not SuSiE negatives. Because the two methods can test different genes, we report each count against the genes that method could test, and we report single-signal-only genes separately. We ran colocalization twice. The second run repeated only the positives from the first run and added the shared-posterior check, which requires both signals in a pair to keep at least half of their posterior probability on variants that the GWAS and eQTL fits share. Negative results from the first run were carried forward unchanged. The count of finite results therefore includes these carried negatives, and we report it separately from the subset that passed the shared-posterior check. Of the 436 primary SuSiE genes, 424 are protein-coding and 12 are lncRNA targets. Gene-level eQTL data do not identify transcript isoforms or give a biotype-specific enrichment denominator.

The COLOC gene counts are not derived from the 360 fine-mapped loci. Of the 436 multi-signal genes, 100 have a colocalized lead variant within 250 kb of a fine-mapped window from the same GWAS (89 within the window itself). Of the same 436 genes, 94 are supported only by the Million Veteran Program (MVP) European ALT, AST, and NAFLD GWAS [ref: MVP liver-trait GWAS], which were not fine-mapped. To describe where the lead colocalized variants fall, we assigned one lead COLOC variant to each of the 996 genes (933 distinct lead variants) and classified it by overlap with transcript annotations (Figure 2C). In this gene-level summary, the lead variant lay in coding sequence or a splice site for 35 genes, in a UTR for 70, in a promoter for 104, in an intron for 524, and in intergenic sequence for 263. Of the 35 coding genes, 14 have the lead variant in their own coding sequence or splice site and 21 only in another gene's. For example, the lead variant of SAMM50 is the PNPLA3 variant, and that of MAU2 is the TM6SF2 variant.

We describe the complete result as a map of inherited variation relevant to MASLD and liver traits, and we reserve direct MASLD language for the tiered disease subset. We restrict regulatory language to supported GWAS-eQTL relationships; protein-altering and unresolved noncoding signals keep their own labels. The eQTL side reuses an existing panel (Broadaway et al., 2024). This analysis therefore adds MASLD and liver-trait tiering, multi-ancestry GWAS integration, fine-mapped DNA context, and comparison with diseased human liver, not a new molecular resource. A colocalized gene is a genetically anchored candidate; colocalization alone does not show that the gene mediates the trait association.

Ancestry labels describe the GWAS populations. They do not mean that molecular regulation was ancestry-matched in every analysis. The liver-eQTL panel is European, so a non-European GWAS paired with it tests cross-population sharing, not ancestry-matched liver regulation. We report LD and eQTL provenance for these pairs and do not treat them as ancestry replications (Figures 2E and 2F). Of the 436 multi-signal genes, 340 were supported only in European studies, 64 in European and at least one non-European study, and 32 only in non-European studies (Figure 2E). Of the 340 European-only genes, 213 were not tested in any non-European study, so their European-only status is untested rather than unsupported. Thirteen genes had PP.H4 > 0.5, by SuSiE or ABF, in four of the five GWAS ancestry groups, and no gene was supported in all five (Figure 2F). These partitions do not establish biological differences between ancestry groups.

Two locus views show what the Resource records for single genes (Figures 2G-2I). At RORA, the MVP European NAFLD GWAS colocalized with the liver eQTL (SuSiE PP.H4 0.983), although no variant in the window reached genome-wide significance (lead *p* = 7.6 x 10^-8) and this GWAS has no fine-mapping record (Figure 2G). At FABP1, the UK Biobank ALT GWAS colocalized with the liver eQTL (SuSiE PP.H4 0.998), but FABP1 was not differentially expressed in bulk liver (log2 fold change −0.04, *p*_adj = 0.51) (Figure 2H).

We next asked which trait class supplies each gene's strongest colocalization evidence. For 423 of the 436 genes it came from a liver-enzyme trait (ALT 137, AST 158, GGT 128), and for 13 from a direct disease or liver-fat trait (NAFLD 8, PDFF 4, NASH 1). Each gene is counted once, under the trait class of its strongest colocalization, so this split does not show how many genes each trait class supports. Counting any direct-trait support, 22 genes have a direct-trait SuSiE PP.H4 above 0.5: 12 direct only, 11 of them SuSiE-tested in an enzyme GWAS, and 10 with both classes. Among the 362 genes SuSiE-tested in both trait classes, 56 are enzyme-supported and 21 direct-supported (Figure 1C). Multi-signal colocalization tested far more enzyme than direct-trait gene-by-GWAS pairs (12,121 against 771), of which 739 (6.1%) and 28 (3.6%) colocalized. These per-pair rates are descriptive. Pairs within a locus are not independent, and we made no test. The second run repeated only first-run positives (816 gene-by-GWAS pairs). The 171 direct and 3,844 enzyme credible-set pairs eligible in that rerun therefore do not measure per-signal yield, and 124 of the 171 direct pairs come from MVP European NAFLD, which has no fine-mapped credible sets. Most colocalized genes thus take their strongest evidence from liver-enzyme traits, a pattern that largely follows how many gene-trait pairs were tested in each trait class.

## The Resource assigns each gene a genetic, disease-state or two-axis evidence class

Having defined the colocalized gene set, we compared it with the disease-state transcriptome from the bulk fit over 774 genetic-identity clusters, which are not verified donors (Figure 1B). Among the 3,769 stable-Ensembl genes evaluated by multi-signal colocalization, 426 were SuSiE-colocalized and 220 were canonical disease-state genes. The two sets shared 26 genes against 21.2 expected within expression deciles (Mantel-Haenszel OR 1.29, 95% CI 0.84 to 2.01; exact *p* = 0.28, one prespecified test). Expression level was not fully controlled, so we state no direction. The count is a description, not evidence that the two layers are biologically separate or convergent. A further 10,928 jointly tested genes, including 935 canonical disease-state genes, have no multi-signal COLOC result. Another 7 are SuSiE-untestable because, in one tier-1/2 study, no signal pair passed the shared-posterior check. Neither group is counted as negative. Another 699 canonical disease-state genes, 503 of them lncRNAs, have no tier-1/2 colocalization record and lie outside the 14,704 genes with both a bulk test and a colocalization record. Together with the 935 disease-state genes above that have no multi-signal result, this leaves 220 of the 1,854 calls in the comparison. The raw categorical universe contains 14,720 genes and the annotation-adjudicated Catalog universe 14,682; their denominators are not interchangeable.

Across the 3,769 genes, PP.H4 and the absolute disease-state log2 fold change had Spearman rho = 0.001 (nominal *p* = 0.95). Adding the seven SuSiE-untestable genes gives the 3,776 genes drawn in Figure 3G, with rho = 0.0014 (nominal *p* = 0.93). Because PP.H4 is unsigned, we did not test it against the signed fold change. Separating inherited genetic evidence from established tissue state has two practical consequences. First, differential expression cannot serve as a required validation filter for genetically anchored regulatory candidates. Second, absence of genetic support does not make a disease-state program unimportant, because common-variant colocalization has limited power and tests inherited regulation rather than tissue activity. We therefore assign each gene to an evidence class rather than treating either map as a noisy approximation of the other. Genetically anchored genes have SuSiE PP.H4 > 0.5 and no disease-versus-control RNA call. Established-state associated genes have an RNA call and no SuSiE PP.H4 > 0.5. Genes supported on both axes meet both thresholds; the direction of the genetic effect is not compared with the RNA direction. Unresolved genes are testable in both domains and meet neither threshold, and untested genes are absent from, or untestable in, at least one domain (Methods).

We next checked the colocalization state of the recurring disease-remodeling examples. None of GNMT, MAT1A, CYP2C19, IGFBP7, and BICC1 has a primary SuSiE PP.H4 above 0.5. Only GNMT has a single-signal ABF PP.H4 above 0.5, in one direct-NAFLD stratum (0.706, European NAFLD GWAS [ref: 2021 European NAFLD GWAS, PMID 34841290]). Its shared-variant anchor is at GRCh37 6:42931261 (rs2296804, SNP.PP.H4 0.259); SNP.PP.H4 is the posterior that a given variant is the shared one. The NAFLD association at this locus is weak: the anchor has GWAS p = 2.6 × 10⁻⁴, and the smallest p within 1 Mb is 1.4 × 10⁻⁴. The colocalization therefore rests mainly on the strong eQTL (p = 2.9 × 10⁻⁶⁶) and on the prior probability that one variant affects both traits (p12 = 5 × 10⁻⁶). In the earlier run of this pair (PP.H4 0.708 at that prior), PP.H4 was 0.33 at p12 = 1 × 10⁻⁶ and 0.83 at 1 × 10⁻⁵. We therefore treat GNMT's genetic support as sub-threshold sensitivity evidence. Multi-signal SuSiE colocalization was attempted for GNMT in that stratum but returned no credible-set result, and its best SuSiE PP.H4 in any primary study is 0.128 (UK Biobank ALT). The strongest stored ABF anchors for MAT1A, CYP2C19, IGFBP7, and BICC1 are at 10:81143387, 10:96292988, 4:57949517, and 10:60583832, respectively, but their gene-trait PP.H4 values remain below 0.5. These coordinates are COLOC shared-variant anchors, not causal variants or GWAS credible-set PIP.

We treat GNMT as a provisional worked locus. At the nonpalindromic variant rs2296805, the G allele lowers liver GNMT expression and aligns with higher direct-NAFLD risk. rs2296805 is a proxy for the rs2296804 anchor (UK Biobank European r2 = 0.998). The anchor is a C/G variant at a frequency near 0.5, so its allele labels do not resolve strand on their own. The nonpalindromic rs9471976 also aligns lower expression with higher risk and has r2 = 0.818 with rs2296805. rs11752813 (r2 = 0.580) is palindromic, and its reported sign is not independently harmonized. Ji et al. (2012) associated rs9471976 and rs11752813 with GNMT protein in 268 liver biopsies and reported allele-dependent reporter activity, which motivates resolving this locus. Unsigned LD and colocalization do not identify a causal variant or show that GNMT mediates the association.

Expression-based colocalization may miss some protein-altering and protein-level loci. Protein-altering variants can act through protein function without changing steady-state RNA, so they need not appear as expression QTLs; others change RNA abundance through nonsense-mediated decay or splicing (Rivas et al., 2015). SuSiE credible sets place a missense variant at three such genes: MTARC1 (rs2642438; normalized posterior inclusion probability, PIP, 0.757 in a European PDFF GWAS), TM6SF2 (rs58542926; PIP 1.0, European PDFF), and PNPLA3 (rs738409; PIP 0.706, MVP Admixed American NAFLD). All three genes are in the unresolved Catalog class, and external protein-QTL evidence for MTARC1 is a source finding (Gobeil et al., 2026b). The Gene Catalog records such members for every gene. In total, 31 genes carry a fine-mapped protein-altering or canonical-splice credible-set member at normalized SuSiE PIP ≥ 0.5 (20 unresolved, 7 genetically anchored, 3 established-state associated, 1 untested, and 0 supported on both axes). This field misses non-canonical splice effects, does not change the eQTL-colocalization class, and does not show that the gene mediates the trait association.

## Noncoding variants carry most fine-mapped posterior probability at liver-trait loci

To measure DNA architecture directly, we used every reliable tier-1/2 SuSiE credible set that fine-mapping produced, not only the credible sets that can be rebuilt from the top SNP of a colocalized gene. Of 6,467 study-by-locus credible-set instances, 6,415 passed the prespecified PIP-sum criterion. Complete GRCh37 reference-oriented VEP annotation [ref: Ensembl VEP] assigned 100% of their normalized PIP to a consequence class. We grouped repeatedly recovered sentinels within 1 Mb and weighted each genomic region equally. Direct MASLD and liver-fat loci then allocated 23.15% of the probability to protein-altering variants, 0.20% to canonical splice variants, 3.16% to synonymous or UTR variants, and 73.50% to other noncoding variants. Liver-enzyme loci allocated 3.87%, 0.00%, 3.60%, and 92.52% (Figure 2D). The grouping limits the weight of repeatedly recovered studies and does not assert LD or statistical independence. When each credible set was instead weighted equally, protein-altering variants carried 15.77% of the posterior for direct traits and 3.26% for liver enzymes.

Fine-mapping coverage differs by stratum. Every lead window of the European direct-trait GWAS that entered fine-mapping has a credible set (47 of 47). The European enzyme arm consists of 82 credible sets from the UK Biobank ALT, AST, and GGT GWAS [ref: UK Biobank liver-enzyme GWAS]. These GWAS were fine-mapped at 6, 9, and 15 loci, which match 8 of their 767 lead windows by locus identifier (23 of 767 within 500 kb by position), and the arm has no credible set at PNPLA3 or TM6SF2. The MVP European ALT, AST, and NAFLD GWAS were not fine-mapped. In the Admixed American (AMR) GWAS, 10 of 37 direct and 84 of 271 enzyme lead windows have a credible set.

Taking protein-altering and canonical-splice variants together, direct and enzyme loci pooled over ancestries carry 23.34% and 3.87% of their posterior probability, respectively (difference 0.1947, 95% physical-cluster bootstrap interval 0.0692 to 0.3399; 233 physical 1-Mb clusters: 19 direct, 228 enzyme, and 14 shared). For this interval, we grouped credible sets into 1-Mb physical clusters and resampled whole clusters 10,000 times, because nearby credible sets are not independent. This pooled value is not a trait-class contrast. Of the direct credible sets, 106 of 207 are European, and 15 of the 19 direct clusters contain a European direct set, whereas 6,126 of 6,208 enzyme credible sets are non-European. Trait class is therefore confounded with GWAS ancestry, LD reference, and fine-mapping coverage.

To separate trait class from ancestry, we stratified the comparison by ancestry (Figure 2D, ancestry panel; 10,000 shared draws over 81 pooled clusters). European studies gave 29.88% against 6.92% (difference 0.2296, 95% interval 0.0550 to 0.4236), and AMR studies gave 6.73% against 4.54% (0.0219, −0.0489 to 0.1065). The prespecified fixed-weight European-plus-AMR estimate (weights 16/26 and 10/26) was 0.1497 (0.0388 to 0.2718). All non-European studies pooled gave 6.99% against 3.68% (0.0331, −0.0383 to 0.1274), a descriptive value with no calibration control. None of these estimates receives a trait-class call. The European estimate failed three checks: coverage (enzyme arm fine-mapped at 8 of 767 lead windows, direct arm at 47 of 47), power (minimum detectable difference 0.304 against the pooled 0.195), and its random-split calibration control. In this control, we split East Asian enzyme clusters (East Asian and African clusters for the European-plus-AMR estimate) at random into two pseudo-classes with the observed cluster counts, 200 times. Both pseudo-classes come from one trait class, so a calibrated 95% interval should cover zero in at least 180 splits, and the 98.333% interval (Bonferroni over the three estimates) in at least 190. The European 98.333% interval covered zero in 188 splits. The AMR and European-plus-AMR estimates also failed this control (141 and 173 splits at 95%; 150 and 184 at 98.333%), and the European-plus-AMR estimate also carries the European coverage failure.

Within the direct class, the three European PDFF GWAS (37 credible sets, 12 clusters) carry 40.14% against 3.87% for all enzyme traits (difference 0.3626, 95% interval 0.1835 to 0.5455). NAFLD plus NASH (170 sets, 13 clusters) carry 12.44% (0.0857, −0.0164 to 0.2096) (Figure S2S). The direct excess is also concentrated. Five of the 19 direct clusters carry 76.0% of the summed direct cluster fractions, and the largest, one PDFF credible set at chr11:823,586, carries 22.5%. The 57 robustness comparisons (study and direct-cluster omissions, five grouping windows, and one LD-block comparison) gave ratios from 4.37 to 8.28 around the pooled reference ratio of 6.03, a descriptive value. These comparisons include no ancestry or trait stratification and are correlated sensitivities, not replications.

The two trait classes in the 23.34% against 3.87% comparison also differ in GWAS size and effect size. Direct-trait credible sets come from smaller GWAS than enzyme credible sets (median effective n 21,893 against 160,000 over all 6,467 credible sets) and carry larger standardized lead-variant effects (median 0.0368 against 0.0183). Supplementary Note 4 reports a sensitivity analysis matched on standardized lead-variant effect, over all 6,467 credible sets including those that failed the PIP-sum criterion. No permutation *p* value is reported. Comparing the trait classes within one ancestry requires SuSiE fine-mapping of all 767 UK Biobank ALT, AST, and GGT lead windows and of the MVP European ALT, AST, and NAFLD GWAS. This comparison does not test whether the pooled coding difference contributes to the predominance of liver-enzyme evidence among colocalized genes, or to which genes appear in both the colocalized and disease-state sets.

Within the noncoding posterior probability, we treated promoter, liver ABC-enhancer [ref: liver ABC enhancer-gene maps], and liver-lineage open-chromatin annotations as overlapping contexts, not as probabilities that must sum to one. After re-annotation by genomic position (Supplementary Note 6), no noncoding posterior probability was left without a context call. Summed over all credible sets of each class, 13.84% of direct-trait and 14.67% of liver-enzyme noncoding posterior probability lies within 2 kb of a TSS, 9.26% and 8.44% in liver-lineage open chromatin, 0.43% and 3.70% in a liver ABC enhancer, and 77.90% and 77.69% in sequence with none of these contexts. These pooled values mix ancestries. AMR is the one ancestry with comparable fine-mapping coverage in both classes. Within AMR, we detected no difference in the share of noncoding probability with no deposited context (direct minus enzyme −0.0957, 99% interval −0.344 to 0.126; Bonferroni over five contexts). Promoter, ABC, and accessible-chromatin differences received no call because their random-split calibration controls failed. Thus, most noncoding fine-mapped probability lies outside the deposited promoter, open-chromatin, and ABC annotations, and the DNA result does not assign a regulatory mechanism to most noncoding signals.

## A bulk-liver accessibility QTL shares the GWAS signal at 191 of 581 colocalized lead signals

We next asked whether some of these signals share a variant with an accessible chromatin region in bulk liver. We took 581 European tier-1/2 GWAS lead signals that colocalize with a liver eQTL. For each, we tested whether a chromatin accessibility QTL measured in bulk liver from 138 donors (Currin et al., 2025) shares the same GWAS signal (PP.H4 > 0.8 and LD r² ≥ 0.5 between the two lead variants). Of these signals, 191 share a signal with a bulk-liver accessible region, and Currin et al. (2025) already reported at least 45 of them for the same trait and accessible region. Another 315 map to a different variant, 60 had no testable accessible region, 9 colocalized without confirmed lead LD, and 6 were indeterminate. Of the 436 primary SuSiE genes, 143 have a colocalized accessible region. At 108 of them, the same accessible region also colocalizes with the gene's eQTL. The other 35 genes are reached only through the GWAS-eQTL pair. Two pairwise colocalizations with one GWAS signal do not show that the accessibility QTL and the eQTL share a variant, so we do not link these 35 genes to the region. The eQTL and accessibility-QTL donor sets may overlap, so the three-way support is not independent replication.

To estimate how often a shared accessible region occurs, we excluded chromosomes 1 and 19, which a pilot analysis had examined. Of the remaining 238 independent liver-trait loci, 80 (0.336, Wilson 95% CI 0.279 to 0.398) carry a shared accessible region. This rate is conditional on prior liver-eQTL colocalization and is not the share of liver GWAS loci with an accessibility QTL. We then compared loci selected the same way that carry signals from one trait arm only (172 liver-trait loci and 78 platelet-count loci). At equal GWAS and accessibility-QTL strength, liver-trait signals shared a liver accessibility QTL more often than platelet-count signals (adjusted odds ratio 3.35, 95% CI 1.62 to 7.40; *p* = 0.0016, one prespecified test). The unadjusted counts, 70 of 172 against 11 of 78, are descriptive. This agrees with the trait-level pattern reported by Currin et al. (2025), and we do not present it as new.

Taken together, these analyses define the genetic arm of the Resource as a bounded map of trait-associated inherited variation, fine-mapped molecular context, and colocalized liver-trait regulatory candidates. This map does not substitute for the established disease-state transcriptome, which we describe next.

## Five-cohort bulk RNA-seq defines the established disease-state layer

To characterize established disease tissue, we fit a pooled disease-versus-control model to fragment counts from five public liver RNA-seq cohorts [ref: source studies for GSE213621, GSE135251, GSE162694, GSE130970 and GSE126848]. The unit of inference is the genetic-identity cluster, with one library per cluster. These clusters are not verified independent donors. The model includes 774 clusters (138 control and 636 disease; 326, 207, 140, 71 and 30 in GSE213621, GSE135251, GSE162694, GSE130970 and GSE126848), tests 24,115 Ensembl gene rows and has fixed cohort and inferred-sex terms (Figures 3A and 3B). It gave 1,854 canonical calls at Benjamini-Hochberg (BH) padj < 0.05 across all 24,115 genes and |log2FC| > 0.5. We keep colocalization and RNA measurements as separate evidence throughout.

To define what a missing call means, we permuted disease labels within cohort and then added known log2 fold changes to the 774-unit counts (50 replicates). The pooled fit called a uniform 1.0 log2FC change in 50.7% to 100.0% of genes per expression decile and called 0.0% of unchanged genes. We call the smallest change detected with 80% recall the detection limit. It was 0.649 log2FC in decile 8, the lowest of any decile, 0.653 in the most highly expressed decile and above 1.5 in the least expressed decile. In the two lowest expression deciles, a missing call therefore does not exclude a 2-fold change. Those genes are below the detection limit, not unchanged.

We next fit each cohort separately with the same 24,115-gene family and canonical threshold (Figure 3C). The union of per-cohort calls contains 13,761 genes: 9,065 called in one cohort and 4,696 in two or more. Only 2,831 of the 4,696 have the same sign in every cohort that calls them (1,933 in two cohorts, 636 in three, 218 in four and 44 in all five). The other 1,865 (39.7%) are called up in one cohort and down in another, and 1,632 of these include a GSE126848 call. Without GSE126848, 346 of 2,565 multi-cohort genes (13.5%) have mixed signs. These 20,054 cohort memberships are a different set from the 1,854 pooled calls, and shared membership is not replication. GSE135251 membership cannot be matched for genetic ancestry, because only 2 of its 10 controls are of European genetic similarity and 6 are of African similarity.

GSE126848 calls track sequencing depth as well as disease. Its model contains 23 control and 7 disease clusters, all male, and gives 6,147 calls. Direction agrees with the pooled fit for 84.36% of the 1,854 pooled calls. The source identity screen removed 24 of 31 disease records, including all ten women, and 2 of 25 controls. This unequal eligibility may limit how representative the cohort is. The retained disease libraries also have about half the fragment count of the controls on the tested genes (median 19.5 versus 35.3 million). Because voom adds 0.5 to every count before computing log-CPM [ref: voom], a gene with no reads has a higher log-CPM in shallower libraries. Of the 6,147 calls, 437 are such genes, with zero reads in all 30 libraries. Adding log2 library size as a covariate keeps 2,025 to 3,056 of the 6,147 calls, depending on whether library size is counted before or after normalization, and keeps none of the 437 (Supplementary Note 7). We therefore do not read GSE126848 per-cohort calls as disease replication.

Direction agreement among the 1,854 pooled calls ranges from 78.3% to 100% across the five cohorts (Figure S3C). When we removed each cohort and refit the full model, 66.7% to 96.8% of calls kept the same sign, and removing GSE213621 changed the most calls. Leave-out fits on a fixed gene set gave all-gene Spearman correlations of 0.870 to 0.995 and threshold retention of 70.7% to 96.6% (Figure S3C). Neither procedure is independent replication of a pooled result. In a spike-in calibration on these counts, an effect placed only in GSE213621, which carries 45.1% of the summed quality weight, was called in the pooled fit in 87.45% of replicates. An effect placed in any other single cohort was called in at most 0.25% of replicates (1.0 log2FC in both cases). Keeping a call in all five leave-one-cohort-out refits therefore does not show that a real gene's call is independent of any single cohort. Adjusting for log library size, top-100-gene share and detected genes together kept 1,563 (84.3%) or 1,569 (84.6%) of the 1,854 calls under two library-size definitions. Across these two fits and four fits that add one library property each, all-gene log2FC Spearman correlations with the primary fit were 0.979 to 0.992 (Supplementary Note 7). Whether these library properties are technical or reflect disease biology is not resolved.

Within the same cohorts, controls include more units of African genetic similarity than disease biopsies do (31 of 138 against 37 of 636). We added four genetic ancestry principal components (PCs), computed from RNA-called genotypes, to the disease-versus-control model. This changed 233 of the 1,854 calls (140 lost, 93 gained), more than within-cohort-and-group permutations of the same PCs (Holm-adjusted *p* = 0.0199 over two tests). Restricting the fit to units called as European genetic similarity (EUR) changed 570 calls (Holm-adjusted *p* = 0.0199) (Supplementary Note 7). The changes are threshold crossings. No lost call changed sign, and 127 of the 140 lost calls fell to |log2FC| at or below 0.5 with padj still below 0.05. Over the 1,854 calls, log2FC with the PCs correlates with the primary values at Pearson 0.9972 (descriptive). The call count is therefore sensitive to ancestry near the 0.5 threshold. In total, 64 genes, GSTM1 among them, met the ancestry-sensitive rule. This label is descriptive, and the primary fit and its 1,854 calls are unchanged.

To test whether disease-state genes also vary with liver fat outside biopsy cohorts, we used GTEx liver [ref: GTEx]. In GTEx, fat was measured on each donor's liver H&E slide as the clear-vacuole area fraction [figure panel to be assigned]. We tested 1,002 disease-state calls measurable in GTEx liver, in 222 postmortem donors not selected for liver disease. Across these genes, changes with image-measured histologic fat went in the same direction as the disease-versus-control changes in biopsies (rho 0.638). This held in the 138 donors who did not die on a ventilator (rho 0.632). The genes also followed fat more than ischemic time (rho 0.636 for fat and -0.350 for ischemic time in a joint model). Two negative controls did not meet their predictions. Slide tissue area gave rho 0.354 (*p* = 0.056; predicted |rho| < 0.2), and ischemic time gave rho -0.350 (predicted -0.2 to 0.3). The prespecified stop rule for tissue area required both *p* < 0.05 and |rho| of at least half the fat rho (0.319). Only the second condition held, so the rule did not stop the analysis. The fat measured here is histologic steatosis of unrecorded cause in postmortem liver. This result does not assign any biopsy-cohort call to steatosis.

We then scored the adult signature in two pediatric cohorts [figure panel to be assigned]. The adult signature score is built from the 1,854 adult disease-versus-control calls, and 1,064 of their gene symbols map to the GSE66676 array. GSE66676 profiled liver from severely obese adolescents at bariatric surgery, so its controls are also severely obese [ref: GSE66676 source study]. In GSE66676, the score did not differ detectably between 32 adolescents with NAFLD and 33 without (standardized difference d 0.44, 95% CI -0.07 to 0.95). The cohort had 80% power for d of 0.71, and the median adult NAFL-versus-control shift (d 0.97) lies outside this interval. In 52 children in GSE185051 [ref: GSE185051 source study], the score correlated with NAS at rho 0.12 (nominal *p* = 0.41).

We next asked whether a prespecified phenotypic endpoint, rather than groups chosen to improve unsupervised separation, reproduces the disease-state comparison (Figure 3D). The endpoint compares 113 strict controls with 279 advanced-disease clusters across four cohorts. Advanced disease means NAS of at least 5 or fibrosis of at least F3. Strict controls have no observed NAS above 1 and no observed fibrosis above F0. Missing histology is allowed only when no observation contradicts source control status. Each cohort must provide at least five units in each arm. Expression is adjusted for cohort and inferred sex without protecting endpoint labels, and the density and scatter plots show the same 392 clusters. Across the complete 24,115-gene family, endpoint effects correlate with pooled effects at rho = 0.953. Under the same threshold, applied independently, the endpoint recovers 1,755 of 1,854 canonical calls and makes 1,195 endpoint-only calls. This is a related sensitivity analysis in overlapping units. It does not replace the five-cohort disease-versus-control estimand.

Multicohort differential expression and stage-ordered transcriptomic summaries in MASLD are established analyses, not contributions of this paper. A nine-cohort public bulk meta-analysis of 897 NAS-scored and 856 fibrosis-scored biopsies reports stage-resolved differential expression and pathway change on overlapping data, including three of the five cohorts used here [ref: Kumak et al., 2026, bioRxiv, doi:10.64898/2026.08.26.747338]. It follows earlier harmonized liver resources [ref: earlier harmonized liver transcriptome resources]. We therefore claim no novelty for recovering a multicohort disease-state gene set or for ordering expression change by histologic stage. What Figure 3 adds is the estimand and the unit: one fragment-count dataset, genetic-identity clusters as the inferential unit in every bulk fit, a complete 24,115-gene family carried through each contrast, and the connection of these changes to the inherited map (Figures 3F and 3G). The five cohorts include no repeat-biopsy cohort. Where the Resource uses repeat biopsies, they enter only within-person paired models and never a cross-sectional contrast.

## Expression differs from F0 at every fibrosis stage, and the F1 difference depends mainly on one cohort

To describe cross-sectional stage-associated remodeling, we compared each fibrosis stage from F1 to F4 with a common F0 reference (Figure 3E). We used the three cohorts that record an exact F0–F4 stage (GSE130970, GSE135251 and GSE162694; 387 staged units, F0/F1/F2/F3/F4 = 103/102/88/69/25). GSE213621 records only F0–F1, F2 and F3–F4 bins and is not used in these contrasts. Each model has fixed cohort and inferred-sex terms when estimable. BH correction covers all 24,115 genes, and there is no stage-specific fold-change filter, so these counts are distinct from the canonical disease-versus-control calls. The contrasts gave 1,338, 4,266, 6,190 and 7,043 genes at BH q < 0.05 for F1, F2, F3 and F4 versus F0. These counts depend on sample size, effect and uncertainty, and they cannot locate where remodeling begins. Figure 3E also shows continuous F0-centered activity for the five available axes of a six-factor non-negative matrix factorization (NMF) of bulk expression (k6 NMF; the sixth axis, k6-P6, is unavailable). These axes describe continuous variation, not patient classes. Full fibrosis, NAS and composition views are in Figure S3 [panel letters to be assigned].

The F1 count depends mainly on GSE130970. Without it, 14 of the 1,338 genes remain (155 units). Fitted alone, the cohorts give 252 (GSE130970, 50 units), 3 (GSE135251, 90 units) and 0 (GSE162694, 65 units) genes (Figure S3 [panel letter to be assigned]). In a cohort-by-stage model, 95.8% of the 1,338 effects have the pooled sign in all three cohorts. The cohorts therefore mostly agree on the direction of the F0-to-F1 differences, but only GSE130970 gives many calls on its own. The F2, F3 and F4 versus F0 counts stay between 1,934 and 9,546 when any one cohort is omitted. The later adjacent steps do depend on single cohorts. F1 to F2 gives 0 genes without GSE162694 (133 units), and F2 to F3 gives 2 without GSE135251 (54 units). Effect sizes also differ between cohorts at F4, where 4,047 of the 7,043 F4-versus-F0 genes have a cohort-by-stage interaction at BH q < 0.05. We prespecified a test of whether within-cohort library size, library Gini coefficient and normalization factor account for the fibrosis- and NAS-stage calls. We could not interpret that test, because its positive control, a planted covariate, failed. The Figure 3E counts are unchanged.

We also asked whether grouped fibrosis-stage effects from these three cohorts appear in GSE213621 [figure panel to be assigned]. This check is not blind. GSE213621 stage labels already enter the continuum analysis. We had also read a fit of the same F3–F4 versus F0–F1 contrast in GSE213621, and its sign agreement with the discovery cohorts, before the analysis plan was written. The result therefore describes the size of the agreement and its interval. It is not a prediction confirmed on unseen data. We refit the three cohorts on the same bins in disease units only (371 units), with controls excluded on both sides. Grouped effects point the same way in the 258 GSE213621 disease units. The slope was b = 0.708 (95% bootstrap interval 0.550 to 0.867) over the 8,449 genes that have BH q < 0.05 for F3–F4 versus F0–F1 in the three-cohort fit and pass the GSE213621 expression filter. This is larger than the slope that submission order and sex structure produce (permutation null mean 0.006; Holm-adjusted *p* = 0.002). Whether this agreement is as large as the agreement among the discovery cohorts is not determined (ratio 0.793, 95% interval 0.488 to 1.098). The F2 versus F0–F1 contrast was untestable, so we claim nothing about F2 (Supplementary Note 7). A negative control failed. The GSE213621 male-versus-female coefficient is also associated with the discovery stage effects (slope -0.155, 95% interval -0.279 to -0.031, which excludes zero). Either sex and fibrosis bin are not separated in GSE213621, or sex and fibrosis share real biology, for example in sex-biased hepatocyte metabolic genes. Figures 3E and 3F remain three-cohort results.

We kept histology as the clinical severity anchor and asked a separate question: whether an externally anchored transcriptomic continuum could resolve molecular variation that fibrosis categories compress (Figures S3I–S3U). All five cohorts were scored and enter descriptive continuum-window displays, and the four cohorts with recorded fibrosis enter stage displays. We evaluated two prespecified co-primary scores inferentially in GSE162694 and GSE213621, with fibrosis stage and inferred sex in every genetic-cluster-level model. The continuum recovered the broad stage-associated direction of remodeling in transcripts, pathways, NMF axes and Hotspot programs (117 prespecified single-cell gene programs, each a weighted set of co-expressed genes; Figure 4) and kept signal after adjustment for stage (Figures S3M–S3U and S4H–S4N). For the Hotspot programs, that signal depends on cell-type composition adjustment. In the 367 disease units, without the 68 GSE213621 controls, 94 of the 116 estimable programs (of 117) are continuum-associated (BH q < 0.05) without adjustment and 53 to 76 after adjustment, depending on the adjustment arm and the deconvolution input. These comparisons use the same expression data, so they are alignment analyses, not independent validation. Sliding windows are fixed descriptive summaries only. This extension does not change Figure 3, the definition of the canonical disease-versus-control calls or the interpretation of NMF axes.

## GNMT, MAT1A and CYP2C19 have lower RNA and protein abundance in disease liver

GNMT, MAT1A and CYP2C19 are three hepatocyte metabolic genes that illustrate the disease-state result. All three are lower in the pooled disease-state comparison and in adjusted liver protein (Figure 3E). The F0 group contains 14 source controls and 89 disease biopsies, so we also refit the stage contrasts with a disease-only F0 and with an added control-status term. All three genes have lower abundance at F2, F3 and F4 than at F0, with BH q < 0.001 under every F0 reference and without GSE130970. Omitting GSE135251 or GSE162694 leaves every effect negative. Every q is ≤ 0.034 except MAT1A at F4 without GSE135251 (q = 0.059; 11 F4 units). The differences at earlier steps depend on single cohorts (Supplementary Note 7), and we do not use them to place where the lower abundance begins. MAT1A falls below the canonical pooled 0.50 effect floor. It passes the interval-null sensitivity, a test that the absolute log2FC exceeds 0.25 (BH q < 0.05). We do not relabel it as a canonical call. Lower bulk abundance of these genes can reflect a smaller hepatocyte share of the biopsy as well as lower per-cell abundance. These data support a recurring metabolic abundance pattern. They do not show a decline over time within a person, a shared upstream regulator or a patient subtype.

Within 46 MASLD participants with liver proteomics [ref: PXD051911 source study], the three proteins also correlate positively with each other after adjustment for acquisition batch, age, BMI and sex (rank-residual partial Spearman rho = 0.406–0.472 over the three pairs; three-pair BH q = 0.00243–0.00772, accounting for nuisance-model degrees of freedom) (Figure S3 [panel letter to be assigned]). All adjustment comparisons use the same 46 participants. Adding NAS reduces the correlations to 0.223–0.277. Paired participant-bootstrap intervals for these reductions exclude zero, but neither the reduction nor a failed significance threshold establishes mediation (Supplementary Note 7). This example is selection-conditioned, because the protein cohort contributed to its selection. It supports a recurring metabolic abundance pattern, and independent liver-protein confirmation is unresolved.

## No association was detected between colocalization support and disease-state RNA change across 3,769 genes

Figure 3F places stage-associated expression beside regulatory genetics and keeps the two measurements separate. Fifty-four genes belong to the 436-gene primary SuSiE set and pass BH FDR < 0.05 in at least one adjacent fibrosis contrast. The adjacent contrasts (F0 to F1, F1 to F2, F2 to F3 and F3 to F4) contain 103/102, 102/88, 88/69 and 69/25 units (205, 190, 157 and 94 in total), and no identity cluster enters both sides of a contrast. Each contrast is its own BH family, so the union carries no joint 5% FDR. When the 5% budget is split across the four contrasts (BH q < 0.0125 in each), 15 genes pass (11 / 2 / 2 / 1 marks). This bounds the FDR over all discoveries in the four 24,115-gene families at 5%, but not specifically within the SuSiE subset. The four complete 24,115-gene families give 216 displayed effects and 59 significant marks (30 / 13 / 11 / 5). Twenty-nine of the 54 genes are selected only through an F0-to-F1 mark, and those marks share the cohort dependence of the F1 count. None of the 30 keeps BH q < 0.05 without GSE130970, and 18 keep it with disease-only F0. Later marks also depend on single cohorts. None of the 13 F1-to-F2 marks keeps q < 0.05 without GSE162694, and 1 of the 59 marks (OMG, F2 to F3) keeps it after each single-cohort omission. Each omission also removes 12% to 69% of a contrast's units, which lowers power. For each gene, the genetics column shows the PP.H4 of one exact gene-trait pair: the gene with its strongest eligible driving GWAS trait. Fine-mapping PIPs are not shown. Gene selection is conditioned on stage significance, so the panel establishes neither enrichment nor causality. Each gene keeps its evidence state, linkage alternative and next discriminating experiment.

Colocalization does not place GNMT, MAT1A or CYP2C19 in an inherited regulatory program. None of the three has a primary SuSiE PP.H4 above 0.5 (GNMT 0.128, from the UK Biobank ALT GWAS). As described for Figure 2, GNMT's only colocalization above 0.5 is a single-signal ABF PP.H4 of 0.706 in a European NAFLD GWAS [ref: 2021 European NAFLD GWAS, PMID 34841290]. Multi-signal SuSiE colocalization returned no result for GNMT in that study, as for 18,939 of the 18,961 gene rows in that study. There is therefore no multi-signal estimate to confirm or contradict the single-signal result. GNMT is a provisional link between the genetic and disease-state maps. The RNA, stage and protein observations stand independently of that result.

Because PP.H4 is unsigned, we compared colocalization support with the magnitude of disease-state RNA change (Figure 3G). Across the 3,769 genes evaluated by multi-signal COLOC (Figure 1B), Spearman rho between PP.H4 and |log2FC| was 0.001 (nominal *p* = 0.95). With the seven SuSiE-untestable genes added (3,776 genes), it was 0.0014 (*p* = 0.93). The 426 genes with PP.H4 > 0.5 split 218 up and 208 down. The panel groups Neither (unresolved), Transcriptomics only (established-state associated), Genetics only (genetically anchored) and Both (supported on both axes) contain 3,149, 194, 400 and 26 genes. Without stratification by expression, Both holds 26 genes against 24.9 expected under independence (Fisher OR 1.06, nominal *p* = 0.83). The prespecified test, stratified by expression decile, gives 21.2 expected (Mantel-Haenszel OR 1.29, *p* = 0.28; Figure 1B). The seven SuSiE-untestable genes are drawn separately and are not counted as Neither. This gene universe is distinct from the annotation-adjudicated universe of the Gene Catalog and is not replaced by it.

We separately compared 426 one-to-one RNA-abundance-matched pairs, each a genetic gene and a control gene, using the same interval-null sensitivity (lfc = 0.25, BH q < 0.05) (Figure S3G). Calls occur only in the genetic member in 37 pairs, only in its control in 28, in neither in 360 and in both in one (paired OR = 1.321; nominal exact McNemar *p* = 0.321). Genes are not biological replicates. Because PP.H4 is unsigned, neither comparison tests risk-allele direction or establishes biological opposition or complementarity.

The inset of Figure 3G compares lncRNA with protein-coding genes in the 774-unit RNA fit. Over 21,501 genes, the odds of a call for lncRNA relative to protein-coding genes were 1.12 without adjustment (95% CI 1.00–1.25) and 0.43 after adjustment for expression and detection (95% CI 0.38–0.49). Calls also differed in direction: 83.0% of 1,157 protein-coding calls and 47.9% of 549 lncRNA calls were up. We matched 520 one-to-one lncRNA/protein-coding pairs on covariates at an outcome-blind 0.15-SD caliper. In these pairs, the conditional odds ratio for upregulation, lncRNA versus protein-coding, was 0.132 (30 discordant pairs with the lncRNA member up against 227 with the protein-coding member up; Wald 95% CI 0.090–0.193). Gene-level matching is descriptive, not biological replication.

Having described the bulk disease-state and stage layers, we next examined prespecified cell-type-specific gene programs and regulatory activities at the level of biological donors (Figure 4).

## A 102-donor single-cell atlas contains 117 prespecified cell-type-specific gene programs

To place disease-associated gene programs in specific liver cell types, we integrated seven single-cell datasets (275 libraries, 1,232,318 cells) into one annotated atlas [ref: source studies of the seven single-cell datasets] (Figure 4A). After two prespecified exclusions of low-cell libraries, the analyzed atlas has 273 libraries (1,232,285 cells) from 102 donors. It carries 117 prespecified cell-type-specific gene programs, weighted sets of co-expressed genes [ref: Hotspot], across five tested lineages. We summarized program, communication and regulator analyses at the donor level, because sequencing runs, sort fractions and cells are not biological replicates. Primary-source diagnoses give 24 Healthy, 14 Steatosis, 36 Steatohepatitis, 17 Cirrhosis and 11 Unresolved donors. The GSE174748 source study identifies its two NAFLD explants as cirrhotic [ref: Filliol et al.]. Table S6 of the GSE192742 source study resolves its donors as 4 Healthy, 3 Steatosis, 1 Cirrhosis and 11 Unresolved [ref: Guilliams et al., 2022, Table S6], and all 19 GSE192742 donors remain excluded from stage inference. A generic NAFLD or fibrosis label alone does not establish a diagnosis of steatosis or steatohepatitis.

The atlas follows the same evidence rules as the rest of the Resource. A program received an evidence state only when donor identity, retained gene weight, cell-state coverage and complete multiple-testing families were all in place. If any one failed, the program was untestable. Untestable is not evidence that the biological program is absent. Each selected program keeps its unresolved alternative and its next discriminating experiment. The integrated UMAP shows 50,000 cells sampled uniformly from the analyzed cells and colored by cell type (Figure 4B). It is descriptive, and cells are not the unit of inference. We do not draw a donor-averaged disease density, because its Healthy and Steatohepatitis groups come from different mixes of datasets.

## Gene content contradicts the discovery lineage for 37 of 117 programs

A program's name records the cell-type subset in which it was discovered, not what its genes are. To test whether the genes match the name, we scored each program against donor-level expression in each lineage. Content matched the declared lineage for 80 of 117 programs and contradicted it for 37 (Figure 4D). Across 1,000 donor bootstrap draws, 76 matches and 34 mismatches kept the same call in at least 95% of draws, and 7 programs were unstable (10 when donors were resampled within each lineage). The per-draw match count had a central 95% range of 77 to 80. The mismatches are systematic. Of the 37, 17 were discovered outside hepatocytes but have hepatocyte content, and several are named for it. Examples are a "Hepatocyte xenobiotic (CYP2C)" program labelled fibroblast and a "Bile-acid transport (ABCB11)" program labelled cholangiocyte. These 17 programs carry a higher measured burden of ambient RNA, the free transcripts captured along with each cell, than correctly labelled non-hepatocyte programs (median 0.432 against 0.061). Hepatocyte transcripts are the dominant ambient species in these data. This association does not show that ambient RNA explains every mismatch. The panels separate three lineages for each program: the lineage its content matches, the lineage it was discovered in and the lineage used for disease testing. The program definitions are unchanged.

In separate donor-paired contrasts within the same atlas, Stromal ECM/IGFBP7 scored highest in fibroblasts (beta = 2.271 over hepatocytes; 53 donors) and Ductular injury/BICC1 scored highest in cholangiocytes (beta = 1.574 over hepatocytes; 51 donors) [figure panel to be assigned]. Both pass the ten-contrast BH and HC3 families [ref: HC3 heteroskedasticity-consistent standard errors]. This is descriptive localization. It does not show where the transcripts came from or that the program acts within that cell type.

## No program shows a supported ordinal stage association in two multi-stage single-cell datasets

Two focal programs, Stromal ECM (IGFBP7) and Ductular injury (BICC1), were discovered in hepatocyte profiles but are among the 37 whose content contradicts that label. To test all 117 programs for association with ordered disease stage, we coded the two GSE174748 explants as Cirrhosis, as documented in the source study, excluded cirrhotic donors under the primary Healthy/Steatosis/Steatohepatitis estimand and used the 60 donors whose datasets could be ambient-corrected (Figure 4C). Only GSE202379 and GSE244832 contain more than one stage, so with dataset as a covariate they alone inform the slope (55 hepatocyte, 38 cholangiocyte, 43 fibroblast, 25 macrophage and 53 T-cell donors). Residual variance differs between these two datasets. We therefore fit a separate variance per dataset by generalized least squares [ref: generalized least squares implementation] and compared each program with 19,999 within-dataset stage permutations. No program passed the 117-program BH and stability rule on the uncorrected scores (smallest q = 0.323). This rule requires BH q < 0.05 over the 117 programs and a program whose top 50 genes are recovered when each dataset is left out in turn (median best-match Jaccard index at least 0.5). Ductular injury had beta = 0.137 (95% CI 0.012 to 0.261; permutation *p* = 0.037, q = 0.323). The probability of detecting an effect the size of its estimate in the 62-donor equal-variance fit (0.402; Supplementary Note 8) was 0.999. Stromal ECM had beta = 0.084 (−0.067 to 0.235; q = 0.578). The detection probability at its estimate in that fit, 0.352, was 0.79, so this analysis cannot exclude an effect of that size. The two dataset slopes differ for both programs (Ductular injury 0.107 and 0.742, Cochran Q *p* = 0.016; Stromal ECM 0.037 and 0.756, Q *p* = 0.0068). Each pooled beta is therefore a precision-weighted average, not one stage effect. A 62-donor equal-variance fit is reported as a supplementary sensitivity (Supplementary Note 8). This null result bounds only a pooled ordinal stage association in these two datasets. It does not show that the programs are absent.

We next asked how much ambient RNA contributes to the two focal programs. In an ambient-RNA correction analysis of at most 60 donors across four datasets, correction lowered Stromal ECM from raw beta = 0.3509 to corrected beta = 0.1279 (paired delta = −0.2230, a 64% reduction; nominal *p* = 0.0064; 117-program delta BH q = 0.4781). Ductular injury moved from raw beta = 0.3846 to corrected beta = 0.4010 (delta = 0.01642; nominal *p* = 0.67; delta BH q = 1) (Supplementary Note 8). In a descriptive refit of the same per-dataset-variance model to the ambient-corrected scores, outside the primary test, the betas are 0.052 for Stromal ECM and 0.143 for Ductular injury (permutation q = 0.559 and 0.270). This refit covers the 106 non-T-cell programs, because T-cell programs have no corrected score, and no program reaches q < 0.05. The attenuation classification applies only to programs selected in the primary fit, and neither program was selected, so it applies to neither. Together with the gene content of both programs, which contradicts their hepatocyte label, these results do not let us assign either program to hepatocytes in the original atlas.

We then tested the two programs in independent donors from GSE256398 [ref: GSE256398 source study], using single-nucleus hepatocyte pseudobulk. The models were score ~ stage + sex + age and score ~ disease + sex + age, with HC3 standard errors and Holm correction over four tests [figure panel to be assigned]. No test passed Holm correction. For Ductular injury, the stage beta was 0.694 (95% CI 0.062 to 1.326; 13 non-cirrhotic donors; Holm *p* = 0.139). Its disease-versus-healthy beta was 1.080 (−0.091 to 2.250; Holm *p* = 0.203) in 17 donors, and the 11 disease donors in that contrast include 4 with MASH cirrhosis. For Stromal ECM, the two betas were 0.300 (−0.709 to 1.309) and 0.176 (−0.840 to 1.192), Holm *p* = 1. With 13 and 17 donors these tests have low power, and no positive control was run.

We next compared steatohepatitis with steatosis within each program's content-matched lineage, using the cells and donors for which ambient correction was available (Figure S4 [panel letter to be assigned]). Stromal ECM in fibroblasts had a positive HC3-based meta estimate (beta = 0.329, 95% CI 0.098–0.560; 26 donors; two-program BH q = 0.0106). Ductular injury in cholangiocytes was unresolved (beta = 0.0848, −0.135–0.304; 26 donors; q = 0.449). These 26-donor contrasts use only GSE202379 and GSE244832, so the GSE174748 and GSE192742 source-diagnosis corrections do not affect them. Neither result excludes a mixture of substates within the lineage. Cohort-specific and coverage sensitivities are in Supplementary Note 8. Comparisons against a healthy reference are untestable at the primary threshold, because GSE202379 has fewer than three healthy donors in these lineages (1 cholangiocyte and 2 fibroblast reference donors at the 50-cell threshold). A separate analysis of activity versus abundance in the native lineage left 55 programs indeterminate. The content-lineage analysis asks a different question and does not repeat that result. Captured-cell fractions are reported for the same donors and are not tissue cell proportions.

A supplementary cell-cell communication screen used 69 biological donors across four datasets [ref: ligand-receptor inference method] (Figure S4D). Of 55,402 native groups, each a ligand-receptor pair in a cell-type pair, 349 were fitted and 55,053 had insufficient counts per stage. The native screen used BH correction within 27 families of cell-type pair and stage term, then Bonferroni correction across nine cell-type-pair families. In the fixed 13-pair display, six pairs were testable, all with corrected q = 1, and seven were untestable. In a separate fixed screen of 26 program tests, 22 were estimable and two Ductular injury associations were supported, each in 58 donors (Figure S4E). These were hepatocyte→endothelial CDH1→PTPRM (beta = 0.4957, BH q = 3.628 × 10⁻⁶, HC3 q = 1.167 × 10⁻¹⁰) and endothelial→hepatocyte ST6GAL1→EGFR (beta = −0.2987, q = 0.006914, HC3 q = 0.04880). Both keep their direction when each of four cohorts is omitted. These are RNA associations within one atlas. They do not show ligand secretion, receptor activation or a signaling mechanism. Of four selected transcription factor regulons, two are testable. RELA and SMAD3 retain 74.3% and 75.8% of their original L1 weight and are untestable (Figure S4O).

## Bulk program scores show no linear fibrosis or NAS association at BH q < 0.05

To test whether program scores track histology in bulk tissue, we scored all 112 eligible weighted programs in the 364 genetic-identity clusters from the three cohorts of the 774-unit dataset that record both fibrosis and NAS (GSE130970, GSE135251 and GSE162694) (Figure 4E). Within each cohort, we regressed each score on sex, fibrosis stage and NAS. We combined the cohort estimates by random-effects meta-analysis with the Knapp-Hartung adjustment [ref: Knapp-Hartung] and applied BH correction over 224 tests on the two axes. No program met BH q < 0.05 on either axis (smallest q 0.079 for fibrosis and 0.159 for NAS), with or without composition adjustment under either set of cell-type fractions. NAS associations are shown in Figure S4 [panel letter to be assigned]. Adding 93 GSE174478 source records, which lie outside the 774-unit dataset and are not identity-resolved, raises the meta-analysis from three to four cohorts and gives 16 fibrosis calls and one NAS call. After composition adjustment, three fibrosis calls remain with BayesPrism cell-type fractions estimated from the fragment counts used here (fragment-count fractions) and one with fractions estimated from an earlier count matrix of the same cohorts (earlier-matrix fractions). All of these calls depend on the untruncated Knapp-Hartung variance. Each has a between-cohort Q below its expected value, and none remains when that variance is floored at the REML value. Figure 4E therefore shows no program association with linear fibrosis category or NAS that holds across these analysis choices. A linear coding of category does not capture nonlinear remodeling or change within a person. A matched-background comparison of unweighted scores does not establish disease specificity (Supplementary Note 9).

## Program scores covary across bulk livers in two cohorts not used to build the program map

We next asked whether program scores covary with one another across bulk livers [figure panel to be assigned]. In the four discovery cohorts, 2,426 of 6,195 tested pairs of the 113 bulk-scorable programs had a positive partial correlation that passed BH correction, after adjustment for sex, fibrosis stage and NAS and after removal of the leading residual component. Of these, 650 pairs join programs of the same lineage and 1,776 join programs of different lineages. This discovery step used 469 source records from GSE130970, GSE135251, GSE162694 and GSE174478, scored on counts that precede the identity-qualified 774-unit dataset. Its unit is therefore not the genetic-identity cluster used in the other bulk analyses. In bulk tissue, covariance of program scores mixes cell-type abundance, biopsy content and program activity, so it does not show that cell types coordinate. A prespecified test in 258 GSE213621 disease units (genetic-identity clusters) used 2,328 pairs rediscovered on 371 identity-qualified discovery units, and it was untestable. Its program-group permutation control keeps covariance within each program and removes it between programs. This control still gave a mean excess over the non-coordinated baseline of 0.025 (95th percentile 0.081; prespecified limits 0.02 and 0.05), so the run stopped before any GSE213621 outcome statistic was computed. An earlier GSE193066 comparison rests on participant identities that cannot be qualified and is not counted as replication (Supplementary Note 10).

We then tested the 2,426 fixed discovery pairs, with a separate prespecified test and null, in two bulk cohorts that were not used to build or test the program map: GSE281797 (91 participants) and GSE225740 (87 participants), one biopsy each [ref: GSE281797 and GSE225740 source studies]. In both cohorts, the fixed pairs had positive histology-adjusted partial correlations more often than pairs placed at random among programs of the same lineages. This held after removing 1 and after removing 3 residual components. In each cohort the largest of the permutation and bootstrap *p* values was 0.000999, the floor set by 1,000 permutations and 1,000 bootstrap draws, and the null was rejected in both cohorts under Holm correction over the two. With 1 component removed, the fraction of pairs with a positive correlation exceeded the random-placement mean by 0.158 in GSE281797 (0.848 against 0.690; participant-bootstrap 95% interval of the excess 0.132 to 0.182) and by 0.312 in GSE225740 (0.834 against 0.522; 0.270 to 0.323) (Supplementary Note 10). Participant overlap with the discovery cohorts was not checked by genotype in either cohort. GSE225740 identities cannot be checked at all, because it has no reads and no deposited sex, age or batch. GSE281797 is therefore the stronger replication, because it has full covariates and verifiable sex. In GSE225740, the result held after we added three estimated composition components, log library size and the top-100-gene count share to the design. The GSE281797 composition estimate failed its prespecified quality rule and was not read, so the GSE281797 result has no composition adjustment. None of these numbers is pooled with any GSE213621 or GSE193066 value. Estimated composition is not measured composition, so these results do not show that the covariance is independent of cell-type abundance.

To test whether the same pairs covary within tissue, we used Visium spots in 32 HMSMA livers [ref: HMSMA spatial transcriptomics source study] and the 44 programs measurable in every array. Each donor contributed one array. Three arrays inferred to come from controls were excluded. The source study lists 7 controls among its 35 individuals without a key to the arrays, so 4 of the 32 analyzed donors are probably controls that cannot be identified. Program pairs with positive bulk covariance co-vary across spots, but coordinated pairs do not exceed other positively covarying pairs (134 against 80 pairs; excess 0.0116 Fisher z, *p* = 0.124; an excess of 0.0280 or more would have been detected with 80% power). Among cross-lineage pairs that share no genes, the 160 coordinated pairs correlated more across spots than the 372 non-coordinated pairs of the same lineage pairs (excess 0.0257 Fisher z; Holm *p* = 0.0002; positive in 31 of 32 livers) (Supplementary Note 10). However, 292 of those 372 non-coordinated pairs have a gene-disjoint bulk correlation of zero or below. This difference therefore tracks the sign of bulk covariance, not coordination beyond the leading residual component. Each 55-micron spot holds several cells, so this is co-location at spot scale. It does not show expression in different cells or contact between cells. HMSMA has no key linking arrays to clinical records, so it carries no stage information. Neither the bulk nor the spot results show coordination between cell types.

## Stromal ECM and Ductular injury scores increase along a continuous molecular axis after stage adjustment

We then asked whether the two focal programs vary along a continuous molecular axis when recorded fibrosis stage is kept in the model. In the 367 disease units of GSE162694 and GSE213621, both Stromal ECM and Ductular injury increase along the fixed external projection in stage- and sex-adjusted models at the genetic-cluster level (meta coefficients 0.760 and 0.707 score SD per continuum SD). The cohort-PC1 score gives concordant effects. These two continuum scores come from closely related methods (within-cohort correlations approximately 0.999), so they are not independent corroboration. After cell-type composition adjustment, the Ductular injury meta coefficient changes little (0.60 to 0.73 across adjustments and fraction sets). Its GSE162694 estimate, however, falls from 0.952 to 0.477 with the fragment-count fractions and to 0.639 with the earlier-matrix fractions. The Stromal ECM meta coefficient falls to 0.203 to 0.479, depending on the deconvolution input. In GSE162694 alone, Stromal ECM is not supported with the fragment-count fractions (0.172, 95% CI −0.032 to 0.377; BH q 0.233) (Figure S4 [panel letter to be assigned]).

Figure 4F shows descriptive fixed-window score profiles for both programs in all five Resource cohorts and has no inferential panel. The stage- and sex-adjusted estimates are reported separately (Figure S4 [panel letter to be assigned]). GSE126848, GSE130970 and GSE135251 overlap the construction of the continuum and are shown apart from GSE162694 and GSE213621, the two cohorts independent of that construction. GSE126848 is omitted from stage-stratified panels, because recorded fibrosis is unavailable. GSE213621 records Control, F0–F1, F2 and F3–F4, which are not exact F0–F4 stages. Exact-stage estimates and a paired sensitivity to coarsened histology are reported separately for GSE162694. The complete 117-program family, transcript comparisons, pathway collections and ten fixed NMF axes are shown in Figures S3M–S3U and S4H–S4N. Their alignment with the fibrosis-based results supports reading the two together. Histology anchors severity, and the continuum describes a shared axis of molecular variation among samples within the same recorded histologic category. This alignment does not establish a temporal order or replace the fixed program definitions.

## Forty-three gene-membership systems summarize continuum remodeling between single genes and cohort-wide axes

To describe remodeling within a recorded stage by groups of features rather than gene by gene, we grouped 13,092 correlated features into a few dozen auditable systems (Figure 4G). Before any stage or continuum outcome was read, gene membership with the continuum signature genes excluded placed 5,439 nodes in one graph component and 43 size-ranked systems. Six nodes in three disconnected two-node components remain listed but are outside the system summaries. The systems extend the focal view to the complete testable pathway and Hotspot layers. They are not independent. The 43 system scores behave like about 14 independent systems, and their first principal component carries 41% of the variance (466-unit descriptive roster). System scores use the 24,115 genes tested in the expression data, and memberships, weights and coordinates are unchanged. We overlaid stage-adjusted continuum coefficients on that geometry, estimated in the 367 stage-complete disease units (109 GSE162694 and 258 GSE213621). The 68 GSE213621 controls have no recorded histology and are excluded. In the donor-first estimate, we fitted each participant's collection-balanced system score within each cohort and combined the two cohorts by fixed-effect meta-analysis. The feature-first summary instead takes the median node coefficient within each system and collection and averages over collections. Across the 43 systems, the donor-first stage and continuum estimates have Spearman rho = 0.67 (bootstrap 95% interval 0.38 to 0.84), and 34 share direction. The feature-first summaries give 0.63 and 32. Broad molecular summaries therefore align across stage and continuum, while cohort-specific window shapes remain visible (Figures S4Q–S4U). These are descriptive summaries over nested gene sets and overlapping windows. They are not system-level significance tests, a temporal ordering or a universal ontology.

Pathway similarity networks partitioned into modules are an established way to summarize MASLD transcriptomes [ref: Kumak et al., 2026, bioRxiv, doi:10.64898/2026.08.26.747338], and we claim no novelty for that construction. What differs here is the order of operations. We built the 43 systems from signature-excluded gene membership before reading any stage or continuum outcome, and the 117 Hotspot programs keep memberships and weights fixed in advance. The geometry therefore could not have been selected by the outcome it is used to describe. We did not use centrality in a curated protein interaction database to rank genes or programs at any point, because degree in such a database reflects annotation depth as well as biology. An equal-cohort pooled view of all 43 systems is shown beside the five cohort-specific window maps (Figure S4U; Supplementary Note 11).

Figure 4H compares two coefficients for each system: the fibrosis-stage coefficient (score SD per cohort SD of the recorded stage code) and the stage- and sex-adjusted continuum coefficient (score SD per continuum SD). Both are estimated donor-first in the 367 disease units, with 95% donor-bootstrap intervals on both axes. For 22 systems, both intervals exclude zero. The two coefficients estimate different quantities. Two further kinds of test were run on 435 units that included the 68 GSE213621 controls as a separate stage category. The maxT test uses permutations within stage-by-sex strata and compares each system's statistic with the largest statistic across all 43 systems in each permutation, which controls the family-wise error rate. The competitive tests ask whether a system is more strongly associated than random feature sets of the same size, drawn either as connected sets from the prespecified 30-nearest-neighbor feature graph (the connected-region null) or uniformly from all features (the more permissive uniform size-matched null). On that roster, no system exceeds the connected-region null on both co-primary continuum scores (the fixed external projection and cohort-PC1). One system exceeds the uniform null: S01, the largest system (311 features), whose representative terms are integrin binding, extracellular matrix organization, contractile myofibroblast and cell adhesion. Eighteen systems pass maxT in both cohorts on both co-primary scores. The two scores give system coefficients correlated at r = 0.999, so this is in effect a two-cohort criterion. In the 367 disease units, the unadjusted fit calls all 18 (BH q < 0.05 within 43 systems), but 14 of the 18 lose the call in at least one of 11 estimable composition-adjusted fits, and 9 are called in 5 or fewer. A maxT test that combines both cohorts in one permutation family gives 33 of 43 on the 435 units. Because the system scores are strongly correlated and none exceeds the connected-region null, we report these counts as parts of one diffuse shared axis, not as 33 or 18 separate findings.

Retrospective transfer checks in GSE174478 and GSE240729 did not show continuum-fibrosis correlations beyond those of matched random signatures (Supplementary Note 12). These checks were computed on an earlier 844-record version of the bulk data and were not rerun on the 774 units. Further screening of public data did not identify a new, eligible independent validation cohort with verified donor linkage.

Taken together, these analyses place cross-sectional stage-associated remodeling in donor-aware cell-type and continuous molecular context. They do not redefine the program registry, and they do not treat histologic groups or the continuum as states that change over time within a person. Colocalization gives no positive genetic call for IGFBP7 or BICC1, and the unresolved ambient-RNA sensitivity prevents assigning either program to hepatocytes from these data alone. We next asked whether the tissue organization of these programs exceeds assay-matched spatial background (Figure 5).

## Liver protein, chromatin and tissue data place genetic and disease-state evidence in each assay's own units

We next asked how far the genetic and disease-state evidence extends beyond bulk
RNA. We placed both in liver protein, accessible chromatin and intact tissue, and
kept each assay in its own units (Figure 5). Each panel uses its own input set.
Model specifications, correction families and evidence-state definitions are in
Methods, and panel construction and denominators are in the legend. For every
panel, the source table records applicability, evidence state, the unresolved
alternative explanation and the next discriminating experiment.

We first compared RNA and protein changes in the adjusted PXD051911 liver
proteome [ref: PXD051911 liver proteome source study], which includes 58
participants (46 with MASLD and 12 controls) (Figure 5B). The display universe
of 10,044 genes was fixed before the bulk fit of 774 genetic-identity
units, and 3,054 of these genes were measured as proteins. Across these 3,054
genes, RNA and protein changes had Spearman *rho* = 0.458 and 66% direction
agreement (exact values in Supplementary Note 13). Of the 211
protein-significant genes, 188 had the same sign as the RNA change and 23 the
opposite sign. Sign agreement depended on the RNA call. It was 170 of 177 when
the RNA change was itself significant, against 18 of 34, a chance level, when it
was not. We did not read a nonsignificant protein change as falsifying an RNA
call. Across all 5,836 measured proteins, RNA-protein Spearman *rho* was 0.346.
The display universe is partly selected on RNA. *rho* was 0.487 in the 2,390
genes selected from RNA-based analyses and 0.305 in the 664 genes selected on
genetic evidence alone. The RNA-based selection included differential-expression
contrasts from the first bulk fit of the same five cohorts (846 samples), which
preceded the 844-record fit used for the transfer checks and the 774-unit fit.

The protein data extend the GNMT–MAT1A–CYP2C19 example from Figure 3. All three
proteins were lower in MASLD than in control liver (BH *q* < 0.05), in the same
direction as pooled RNA (Figure 5B). The protein model compared 46 MASLD with 12
control participants and had no severity or hepatocyte-fraction term, so it
measures lower bulk protein in disease. In a post hoc analysis that was not
prespecified, all three proteins were lower at higher NAS across the 46 MASLD
participants (partial Spearman *rho* -0.43 to -0.60; BH *q* ≤ 0.007 over 3
proteins and 2 histology scores). Lower bulk abundance can reflect a smaller
hepatocyte share of the biopsy as well as lower abundance per cell. These data
are not evidence that the three genes share one inherited regulatory mechanism.

We next summarized a fixed 25-protein process view, which carries nine canonical
RNA calls and 23 protein calls (Figure 5C). The same protein cohort contributed
to row selection, so this view is descriptive, not independent validation. CD74
is among the nine RNA calls. HLA-DQA1 has a significant RNA change (BH
*p*_adj 0.025), but its log2FC of 0.38 is below the 0.5 threshold for a canonical
call. Its protein change is not significant (log2FC 1.39, BH *p*_adj 0.080). It
is therefore neither an RNA call nor a protein call in Figure 5C.

To ask whether colocalized variants lie in accessible chromatin of liver cell
lineages, we intersected the colocalization posterior of 436 genes with open
chromatin from two liver snATAC sources [ref: two liver snATAC source studies]
across five shared lineages (Figure 5D). These genes form 767 gene-study pairs:
22 genes have a direct MASLD or liver-fat pair and 424 a liver-enzyme pair, and
ten have both. For 12 of 22 direct-trait genes and 167 of 424 liver-enzyme genes,
less than 5% of the posterior lay in peaks in every lineage. We compared the
posterior mass in peaks with the expectation from the tested variants of the same
window. Averaged over liver-enzyme genes, mass in peaks was 2.5 to 3.1 times this
expectation in each lineage. In no lineage, however, did more than half of these
genes exceed their own expectation (94 to 212 of 424). For the 22 direct-trait
genes, the 95% interval of this ratio included 1. Peak overlap does not assign a
variant to a neighboring gene. Absent overlap is unresolved context, not a
negative test.

We then placed all 117 fixed Hotspot programs in chromatin and protein context
without selecting them by disease association (Figure 5E). In total, 31 programs
contained at least one gene from a colocalized pair that placed at least half of
its posterior in open chromatin of the program's own lineage in at least one
snATAC source (genes counted once per primary-assembly gene). In 3 of these
programs (hepatocyte Ion/membrane potential, cholangiocyte FXR/bile-acid and
fibroblast Elastic-fiber ECM), that gene was RORA, which is open through its
direct MASLD pair. In the other 28, genes were open only through liver-enzyme
pairs. RORA is also open in the T/NK Memory/MAIT program, but only through
liver-enzyme pairs. Its direct pair places 0.37 of its posterior in T/NK peaks of
one source, below the 0.5 threshold. In liver protein from 46 MASLD
participants, 22 of 95 scorable programs were robustly associated with a
histologic-burden axis across participants (BH over 95 programs; members counted
once per primary-assembly gene). Of these, 14 were negative, including
oxidative-phosphorylation, xenobiotic and amino-acid degradation programs, and 8
were positive, including fibroblast ECM proteoglycan and TGFB/collagen programs.
This protein cohort did not contribute to program discovery. Neither spatial
program in Figure 5F contains a colocalized gene or can be scored in protein.

Finally, we examined tissue organization for two programs, ECM/IGFBP7 and
ductular-injury/BICC1, labeled Stromal ECM and Ductular injury in the panel
(Figure 5F). Both were chosen in advance from an earlier stage analysis of the
bulk data. Neither passes the 117-program BH and
stability rule for a pooled stage association in the single-cell stage analysis (Figure 4C). In two
spatial transcriptomics sources (Guilliams et al., 2022; Vu et al., 2025), the
prespecified design adjusted spot scores for cell2location hepatocyte abundance,
the discovery lineage of both programs [ref: cell2location].
We then compared residual spatial autocorrelation (Moran's I [ref: Moran's I])
with a matched null built from genes of equivalent expression, detection,
biotype and correlation with hepatocyte abundance. ECM/IGFBP7 exceeded this null
in both sources (excess 0.139 and 0.147). Ductular-injury/BICC1 did not (0.026
and 0.019) and is indeterminate, not negative. A model added afterwards as a
supplementary test of composition dependence, which adjusts for all 16
cell2location abundances and is plotted in Figure 5F, gave the same result
(ECM/IGFBP7 0.040 and 0.050; ductular-injury/BICC1 0.015 and 0.003). Under that
model, residual Moran's I was positive
for ductular-injury/BICC1 in 4 of 4 donors and 9 of 10 arrays. The raw statistic
therefore does not separate the two programs, and no prespecified test compared
them. The matched-gene null holds donors fixed, so it is not a donor-level test.
At the donor level, ECM/IGFBP7 excess was positive in 3 of 4 donors (exact
one-sided sign-flip *p* = 0.19). With four donors, no sign-flip test can reach *p* < 0.05. The four Guilliams et al. donors were sampled during surgery for other conditions: two gastric bypasses for obesity, a cholecystectomy, and a resection of colorectal-cancer liver metastases in a donor who reported six standard alcohol units daily (Guilliams et al., 2022, Table S6). Recorded steatosis ranged from none to 70%, and fibrosis from none to clear pericellular fibrosis. These donors are not a MASLD cohort, and two of them (H37 and H38) also contribute single-nucleus profiles to the atlas, so the spatial and single-cell sources share participants. Methods describe the donor joins, array-level source dependence and
matched-null construction.

Taken together, liver protein shows lower bulk GNMT, MAT1A and CYP2C19 abundance
in MASLD but does not resolve whether this reflects hepatocyte content or change
per cell. Lineage chromatin places colocalized posterior in accessible regions
for a subset of genes without assigning target genes. Of the two spatial
programs, only ECM/IGFBP7 exceeds matched-gene background. Disease association,
cell-of-origin assignment and tissue architecture are separate claims, and each
needs its own calibration. We next organized these results into Gene Catalog
entries, each with a follow-up experiment that could distinguish the remaining
explanations (Figure 6).

## The MASLD Gene Catalog records evidence class, observability and a follow-up experiment for each gene

To make the Resource usable gene by gene rather than as a ranked list, we built
the MASLD Gene Catalog (Figure 6). For every eligible gene and program, a Catalog
entry records the primary gene class, phenotype provenance, assay and biological
unit, testing universe, source dependence and direction where interpretable. It
also records a machine-readable evidence state and the reason for it, the
unresolved alternative explanation and a deterministic next-experiment rule. The
Catalog covers 28,233 gene identities and keeps the eQTL-colocalization call
separate from the bulk RNA call (Figure 6B). The genetic call is supported when
SuSiE-COLOC PP.H4 exceeds 0.5 in any of 35 Tier-1/2 GWAS. The RNA call is
supported at BH *q* < 0.05 and |log2FC| > 0.5 in one limma-voom fit of 774
genetic-identity clusters, which are not verified independent donors. The state
vocabulary distinguishes supported, discordant, indeterminate, untestable, not
applicable, source-dependent and prespecified adequate-negative states. The
Catalog gene calls use three of them (supported, indeterminate and untestable).
No nonsignificant result is converted into a negative call.

A blank evidence cell can mean that the gene was not detected, the source was
underpowered, the assay cannot observe the proposed mechanism, the dataset lacks
an authoritative biological-unit join, or the experiment has not been done
(Figure 6B). These states have different scientific consequences, and the
Catalog records them with the following limits. None of the Catalog's 56,466 gene-level evidence entries is classed as not
applicable. Every entry also carries a source-dependence flag, separate from its evidence
state, so that flag does not separate one gene from another.
The 9,312 genes with no Tier-1/2 COLOC record carry an underpowered-source label.
No power assessment was made for them, so the label marks a missing record, not a
measured lack of power. The Catalog defines no negative-call rule, so no gene is
recorded as tested negative. The genetics call also does not separate evidence
for distinct causal variants from absence of evidence. Among the 816 Tier-1/2
gene-study pairs rerun with a shared-posterior eligibility check, 456 genes have
at least one eligible SuSiE signal pair. Of these, 20 have PP.H3 above 0.8 and
best PP.H4 at most 0.5, the pattern expected when the trait and expression
signals have distinct causal variants. Eighteen of these 20 are called
indeterminate, the same label given to genes with no evidence for either a shared
or a distinct causal variant. The same pattern is common when every Tier-1/2
SuSiE result is counted, including negative results carried forward from the
first colocalization run without the shared-posterior check. Of the 3,918 genes
with any finite Tier-1/2 SuSiE result, a set that does not require a bulk RNA
test and so is larger than the 3,769-gene comparison universe, 3,335 have PP.H3
above 0.8 in at least one stratum and best PP.H4 at most 0.5. Of these 3,335,
3,319 are called indeterminate.

The Catalog also carries annotation fields that describe each gene without
changing its class or route. For each of the 28,233 genes, it records whether
the protein was detected in liver (58 participants) and plasma (143
participants) in the PXD051911 proteomes. It also records whether the gene is
targeted by the CosMx panel of GSE312698 (Boesch et al., 2026), and in which of
nine liver cell lineages its RNA is detected across donors. Where an assay did
not observe a gene, absence is recorded as untestable, not negative. A liver
protein row was quantified for 173 of 1,129 established-state genes, 203 of 400
genetically anchored genes and 11 of 26 genes supported on both axes. Each gene
also carries its position along the central-portal axis in five non-steatotic
live-donor livers (Visium; Yakubovsky et al., 2026) and its hepatocyte protein
zonation (Weiss et al., 2026). Across 1,337 genes measured in both, the two
agreed in direction (Spearman *rho* 0.39, 95% interval 0.34 to 0.44; positive in
each of five donors). Weiss et al. reported a comparison of RNA and protein
zonation against single-nucleus RNA (their Supplementary Fig. 9). Each
disease-state call also records whether the gene changed with image-measured
histologic fat in 222 GTEx postmortem livers. Of the 1,002 measurable calls, 326
did so at BH *q* < 0.05, and 323 of these changed in the direction seen between
disease and control biopsies. This field concerns histologic steatosis of
unrecorded cause in postmortem liver. It does not assign any biopsy call to
steatosis.

On the genetics side, 31 Catalog genes carry a fine-mapped protein-altering or
canonical-splice credible-set member at normalized SuSiE PIP of 0.5 or above,
and 20 of them are in the unresolved class. This field does not change the
eQTL-colocalization class and does not identify the gene through which a variant
acts. Of 581 European liver-trait GWAS signals that colocalize with a liver
eQTL, 191 share a signal with an accessible region in a 138-donor bulk-liver
chromatin-accessibility QTL (caQTL) map (Currin et al., 2025). At least 45 of
these 191 were already reported by the source. At 108 Catalog genes, the same
region also colocalizes with the gene's eQTL. These are annotation fields on
gene records. They create no regulatory-element record and no
credible-set-to-gene link. Constraint, essentiality and secretion fields are
given for follow-up design and do not enter any class.

We illustrate the record with four worked entries: GLP1R, THRB, HKDC1 and MTARC1
(Figure 6C). GLP1R is untested. Its RNA was not detected and it has no Tier-1/2
SuSiE result, so it routes first to the missing decisive assay. THRB has genetic
support (SuSiE PP.H4 > 0.999, UK Biobank GGT) and an indeterminate RNA call,
because its change (log2FC -0.20) is below the 0.5 floor. HKDC1 is supported on
both axes (SuSiE PP.H4 0.988, UK Biobank AST; log2FC +1.07, BH *q* = 8.0 ×
10^-7). MTARC1 remains unresolved on these expression-based axes. The printed
protein-coding routes propose context and function measurements. They do not
convert genetic support into causality, and they do not rank the entries.

The recurring examples of the Resource show further follow-up experiments. GNMT
and CYP2C19 are lower in disease at RNA and protein levels and are
established-state entries, which the Catalog routes to context-matched
perturbation and rescue. MAT1A is also lower at both levels, but its pooled RNA
effect (log2FC -0.485) is below the 0.50 floor of a canonical call. Its Catalog
class is therefore unresolved, not established-state. COLOC retains GNMT's single-signal direct-NAFLD ABF result (PP.H4 0.706), which rests on a sub-threshold GWAS association (p = 2.6 × 10⁻⁴) and does not extend to MAT1A or CYP2C19. GNMT has no multi-signal SuSiE call (best PP.H4 0.128), so its
Catalog genetics call is indeterminate and it takes the same perturbation and
rescue route. Its local regulatory haplotype is unresolved. The ECM/IGFBP7 and
ductular-injury/BICC1 entries separate three routes: ambient-resistant lineage
measurement, matched-null tissue organization and functional perturbation. With
dataset-specific variances modelled, neither program passes the 117-program BH
and stability rule for a pooled stage association. This analysis cannot exclude
the ECM/IGFBP7 effect size estimated in the earlier stage analysis. The original
single-cell atlas lacks ambient-RNA correction. An ambient-RNA sensitivity
analysis corrected five datasets. Two datasets (GSE189600 and GSE192742) could
not be corrected, and 11 T-cell programs were untestable. It does not establish cell of origin.
ECM/IGFBP7 is spatially supported, whereas ductular-injury/BICC1 is spatially
indeterminate. Ductular-injury/BICC1 therefore routes first to a better-powered
donor-resolved tissue test, and both programs need an independent cell-of-origin
assay. PNPLA3 and MTARC1 remain boundary controls for coding or protein-level
mechanisms that expression-QTL analysis need not recover. Each carries a
fine-mapped missense variant (maximum PIP 0.706 for PNPLA3 and 0.757 for
MTARC1), and this leaves its eQTL-colocalization class unchanged. These examples
show how routing works. They are not a ranking of the strongest genes.

The routing layer assigns the smallest experiment that could distinguish the
remaining explanations (Figure 6D). Of the 28,233 genes, 13,551 are untested,
13,127 unresolved, 1,129 established-state associated, 400 genetically anchored
and 26 supported on both axes. Genetically anchored regulatory candidates route
to allele-aware editing with target-expression and locus-linked phenotypes.
Established-state candidates route to context-matched perturbation and rescue.
Candidates supported on both axes need both relationships connected in the same
model. Unobservable candidates first need the missing decisive assay. Transcript
depletion is only one possible perturbation. Increasing gene activity, testing
regulatory-element changes, measuring protein function, lineage restriction and
context-specific challenge remain separate routes. For noncoding nominations,
the route starts from the molecular object. A disease-associated lncRNA routes
to RNA depletion with measurement of the neighboring gene. A colocalized or
ambiguous lncRNA locus needs paired RNA and locus perturbations to distinguish
DNA-element, transcriptional and mature-RNA effects. The Catalog records genes, not
regulatory DNA elements, and does not link credible sets to genes. No gene is
therefore assigned to allele-aware testing of a regulatory element, although
that route is defined.

The Catalog is deposited as machine-readable tables (Data availability;
[accession to be supplied]). Each displayed entry links to its source graph, its
complete tested universe, a data dictionary and its machine-readable
next-experiment rule. A reader can therefore see, for any gene, which layers
were tested, which were untestable and which experiment is proposed.

Across Figures 1 to 7, the Resource records inherited genetic evidence and
disease-state remodeling for each gene, and places both in cell and tissue
context where an assay can observe them. The Gene Catalog is the reusable
output. For each nominated molecular object, it records the biology, the limit
of the evidence and the next experiment.

## Model scores are reported against task-native baselines and each model's training exposure

For a Resource of this kind, a practical question is how much of the mapped
biology a predictive model can recover, and how much of an apparent recovery
reflects data the model saw during its own training. To separate the two, we
scored each model family against its task's own baseline and reported every
score beside a per-cohort record of that model's training exposure (Figure 7A).
We resolved exposure for each study from the model's own audit and development
crosswalk, not from a publication. Where the record was unresolved, we reported
the score as uninterpretable, not as evidence of skill. Of the six
encoders in Figure 7C, Geneformer V2 [ref: Geneformer V2 model paper] had no qualified
clean held-out study, so it has no clean-only comparison.

We first examined broad cell-class mapping in seven held-out single-cell
studies. For the exact TranscriptFormer TF-Sapiens checkpoint [ref:
TranscriptFormer model paper], the complete 644-dataset source inventory
identified five clean studies and two training-seen studies. With the two-layer
head, its macro F1 across all seven studies was 0.9389, a value that includes
both training-seen studies (Figure 7C). On its five clean studies, the same head
scored 0.963 (Figure 7C). Figures 7C and 7F draw each encoder's score and gain
twice, once on all seven held-out studies and once on that encoder's clean
studies. This task maps broad cell classes, not disease stages.

The study-level deltas are means over the held-out studies, so small studies
carry a large share of each mean. One held-out study contributes 2 donors and
149 cells, 0.3% of the evaluation set, but receives one seventh of the weight.
When we weighted studies by donors instead, two of twelve model-by-head blocks
changed the sign of their delta (Figure 7F). We therefore report both means and
make no sign claim that rests on the unweighted mean.

We next tested pretrained sequence models on two variant-level tasks (Figure
7B). The first was a synthetic reporter assay in two HepG2 contexts [ref:
GSE281364 reporter assay source study]. Of seven pretrained sequence-model rows
(six models), one, HyenaDNA [ref: HyenaDNA model paper], exceeded a ridge on
allele identity alone after Bonferroni correction across the seven rows. Its
gain was 0.082 Spearman (95% interval 0.027 to 0.139; Bonferroni-corrected
lower bound 0.009) (Figure 7B). We did not run sequence-based reporter
baselines, so this gain is over a sequence-free control. The second task was
the strength (|z|) of endogenous variant-to-gene associations. Adding either
pretrained embedding to distance lowered the held-out Spearman correlation below
that of distance alone (HyenaDNA -0.022, Nucleotide Transformer [ref:
Nucleotide Transformer model paper] -0.029; intervals exclude zero). Permuting
which embedding went with which variant did not lower the score (Figure 7B).
The labels in this task are marginal associations that were not fine-mapped.
The task therefore cannot separate a model failure from a target that is not
sequence-determined at the variants tested. Taken together, the evidence that
pretrained sequence models add skill on these tasks rests on one reporter-assay
row, measured against a sequence-free control.

No pretrained sequence model was scored on tissue-level severity. We did score
BulkFormer [ref: BulkFormer model paper], a biologically pretrained expression
model, on severity, but all 109 GSE268273 evaluation samples are in its
pretraining set. To evaluate severity, we used the released linear instrument.
It was fitted to 521 sample records from five cohorts (GSE130970, GSE135251,
GSE162694, GSE174478 and GSE240729), and its expression
transform was chosen on GSE268273 [ref: GSE268273 source study]. This fit is
separate from the 774-unit differential-expression fit used for the Resource
bulk analysis (Supplementary Note 14). Only the first three of these cohorts are
in the Resource bulk fit. The model's training pool records no donor key, so
its unit is the sample record, not the genetic-identity cluster of the Resource
bulk fit. In the inferential units drawn in Figure
7D, the Spearman correlation with fibrosis stage, as each cohort codes it, was
0.6631 in 109 GSE268273 participant records (source-declared identity), 0.6680 in 263
GSE213621 genetic identities [ref: GSE213621 source study] and 0.7860 in 81
GSE276114 participants with MASLD [ref: GSE276114 source study]. GSE268273 is
scored on Kleiner stage. GSE213621 and GSE276114 are scored on their own
three-level fibrosis codings. The GSE213621
count excludes 11 identities with conflicting stage labels. In 82 GSE276114
participants with chronic viral hepatitis, the correlation was 0.5443 (Figure
7D). On the published participant records, 299 in GSE213621 and 177 in
GSE276114, the raw Spearman correlations were 0.6382 and 0.6739. The 367-record
GSE213621 analysis that includes controls is a sensitivity analysis in the same
cohort, not another external cohort.

We then compared the instrument with a ridge fitted inside each cohort (native
ridge). In the two additional cohorts, GSE213621 (299 records) and GSE276114
(177 records), the recorded paired differences, instrument minus native ridge,
were -0.040 (95% CI -0.099 to +0.017) and -0.011 (-0.079 to +0.054). These
native fits use expression scaling from the evaluation cohort, so they are
transductive. Their conditional intervals do not establish equivalence. In the
reference cohort, GSE268273, the paired difference was +0.054 (recorded 95%
interval -0.025 to +0.137; minimum detectable effect 0.115). These comparisons
do not establish that any model class is generally superior. We also tested a
nonlinear model with no biological pretraining on the same folds and the same
30 or 100 principal components. Across four comparisons in 109 GSE268273
participants, it differed from ridge by -0.014 to +0.001, and none of the four
differences exceeded the minimum detectable effect of 0.067. An increment below
a minimum detectable effect is not evidence of no effect, so we report it as a
bound, not as a negative result. Across Figure 7, each statistic is placed
against its own null and threshold, and an underpowered arm is labelled
untestable, not negative (Figure 7E).

[Figure 7 model results to be added]

We next asked how much the recorded histology limits interpretation. The
current data do not calibrate its measurement error. Within-participant stage
SDs from source-reported repeat biopsies were 0.4655 in PXD051911 [ref:
PXD051911 source study] and 0.6433 in the 58 source-declared GSE193066 repeat-biopsy pairs (every participant with a second biopsy) [ref:
GSE193066 source study]. In GSE193066, 28 of the 58 pairs agreed exactly.
GSE193066 identities remain unqualified, so the 0.6433 estimate and its score
ICC cannot establish primary within-person reliability or label-error
calibration. The 58-pair score ICC also requires scoring outside the
instrument's released acceptance bounds, because only three pairs have both
biopsies inside those bounds. We therefore treat these values as inputs to
sensitivity scenarios, not as an empirically qualified error range or a
calibrated ceiling. Rank agreement is not calibration to absolute stage. Serial
biopsies can also reflect true change, sampling and shared technical effects.
Interpretation also depends on what each model already saw.

The next discriminating experiment is a within-person one. Cross-sectional
transfer cannot distinguish an instrument that tracks severity from one that
tracks between-participant differences correlated with severity. Paired
biopsies from the same participant can make this distinction. We specify that
test in advance against the measurement bound above, rather than reading the
result off after the fact.

# Discussion

In this work, we present a human MASLD Resource that records, for each gene,
liver-trait genetic evidence and cross-sectional stage-associated disease-state
remodeling, and keeps the assay boundaries of the two separate. The bulk
disease-state layer comes from 774 genetic-identity units in five cohorts. These
units are conservative clusters, not verified independent donors. The genetic
layer comes from eQTL colocalization across liver-trait GWAS. Among the 3,769
genes that had stable Ensembl IDs, were tested in bulk RNA and had a
multi-signal colocalization result, 426 were genetically supported and 220 were canonical
disease-state calls, and 26 were in both, against 21.21 expected at the same
expression level. We detected neither enrichment nor depletion. Because the
matching did not fully control expression, we state no direction
(Supplementary Note 3). A further 935 disease-state calls had no multi-signal
colocalization result, and these are not genetic negatives. The size of this
overlap does not establish biological separation, equivalence or a privileged
intersection of the two maps. It is also not evidence that the underlying
processes are independent. The signed comparison is bounded too. Genome-wide
rank agreement between ALT TWAS z and the disease contrast lies inside ±0.03
without composition adjustment, but not with hepatocyte fraction in the design,
and we state no direction (Supplementary Note 5).

Most colocalized genes take their strongest evidence from liver-enzyme traits,
largely because far more gene-trait pairs were tested in that class. The
strongest native colocalization evidence comes from liver-enzyme traits for 423
of 436 genes and from direct disease or liver-fat traits for 13. Only 22 genes
have any direct-trait colocalization, and far more enzyme than direct-trait
gene-by-GWAS pairs were tested (12,121 against 771). Pooled across ancestries,
protein-altering plus canonical-splice variants carry 23.34% of fine-mapped
posterior probability for direct MASLD and liver-fat traits and 3.87% for enzyme
traits. This contrast is confounded with GWAS ancestry and fine-mapping
coverage, and no ancestry stratum supports a calibrated trait-class call
(Supplementary Note 4). A comparison of the trait classes within one ancestry
needs fine-mapping of all UKBB enzyme lead windows and of the MVP European ALT,
AST and NAFLD GWAS. Some protein-altering variants may act through protein
function without changing steady-state RNA. Others change RNA abundance through
nonsense-mediated decay or splicing (Rivas et al., 2015). One hypothesis for
which genes appear in both sets is therefore that expression-based
colocalization misses risk variants that act through protein function rather
than RNA abundance. The present comparison does not test it, because it does not
hold gene eligibility, phenotype, study power and expression observability
fixed. These results set one requirement for gene lists built this way in MASLD:
the trait that supplies each gene's evidence should be stated. A gene nominated
by ALT is not the same claim as a gene nominated by liver fat.

Matching on expression did not recover the prespecified prediction that
genetically supported genes are more hepatocyte-specific than disease-state
genes. In 281 expression-matched gene pairs, the observed difference was in the
opposite direction (Supplementary Note 3). Every test of that opposite direction
is post hoc, and none is calibrated. The difference comes mainly from the 56
pairs whose disease-state partner was itself tested by colocalization. We
therefore read it as a descriptive difference conditioned on selection, not as
confirmatory class separation or independent biological replication. The
protein-response comparison of 73 pairs gives no support for a class difference,
and it does not establish equivalence. Expression remains a major observability
difference between the classes. Among genes with a single-cell hepatocyte-specificity
value (399 of the 400 genetic-only and 1,124 of the 1,129 state-only genes),
mean bulk expression of genetic-only genes exceeds that of state-only genes by
3.80 log2 units, and 281 of the 399 genetic-only genes obtain an
expression-matched partner. Unequal observability could affect which genes
appear in both sets, but the present comparison does not isolate that
explanation. Germline regulatory effects may act before biopsy, in restricted
cell states or through small persistent shifts. Disease-state expression can
reflect injury, adaptation or composition. Neither class is automatically causal
or irrelevant.

The recurring examples make this distinction concrete. GNMT has a single-signal
direct-NAFLD colocalization (ABF PP.H4 0.706) but no multi-signal call. With
MAT1A and CYP2C19, it is lower in pooled RNA, in cross-sectional fibrosis
contrasts and in adjusted liver protein. In an exploratory analysis of 46 MASLD
participants, the three proteins covary. The covariation remains after
adjustment for fibrosis stage but not after adjustment for NAS or for all four
histology components (BH over 12 pair-by-adjustment tests). This pattern is
compatible with a shared association with histologic severity rather than a
three-gene circuit, but adjustment does not establish mediation. The ECM/IGFBP7 and ductular-injury/BICC1
programs show a different limit. The original atlas lacks ambient-RNA
correction, the ambient-RNA sensitivity analysis cannot resolve all datasets or
programs, and only ECM/IGFBP7 exceeds matched spatial background. Cross-assay
evidence therefore strengthens specific claims, but disease association, tissue
architecture, lineage origin and inherited mechanism remain distinct.

Each component of this Resource has precedent. Staged bulk liver transcriptomes
have been published (Hoang et al., 2019; Govaere et al., 2020; Pantano et al.,
2021), meta-analyzed with gene-level GWAS signal (Piras and DiStefano, 2024),
organized into network modules (Kumak et al., 2026 [ref: bioRxiv
doi:10.64898/2026.08.26.747338]) and ordered along a molecular continuum
(Kamzolas et al., 2026 [ref: doi:10.1038/s42255-026-01543-7]). Single-cell and
spatial atlases resolve disease cell states (Ramachandran et al., 2019;
Guilliams et al., 2022; Boesch et al., 2026). A single-cell multiome assigns MASH
risk loci to liver cell types (Elison et al., 2025), and large GWAS have combined
expression and coding evidence to nominate genes (Sveinbjornsson et al., 2022).
This Resource adds a per-gene Catalog that keeps colocalization results beside
disease-state calls and records whether each layer was untested or tested
without support. Our reanalysis of a 138-donor liver caQTL map (Currin et al.,
2025) places a colocalized accessible region at 191 of 581 European liver-trait
GWAS signals that colocalized with a liver eQTL. At least 45 of these were
already reported by the source.

The same separation of evidence roles has been reported at the protein layer. A
liver proteogenomic study of 504 biopsies found that differential protein
abundance and protein-QTL colocalization nominate largely different gene sets
(Gobeil et al., 2026b [ref: medRxiv doi:10.64898/2026.06.04.26354903]). Of the
five genes that study nominated, our map places GCKR on the primary
expression-colocalization side, MTTP colocalizes only in the
single-signal (ABF) sensitivity analysis (best PP.H4 0.82), ERLIN1 does not
reach that set (best ABF PP.H4 0.47; Catalog class unresolved) and LGALS1 in the established-state class. MTARC1 lies outside our
expression map but is resolved by the protein-QTL analysis of that study. GCKR
has differing expression-QTL and protein-QTL directions. These cases show the
difference between protein and expression observability. They do not establish
a binary split into causal and reactive genes. A single-cell multiomic liver
study likewise reports that chromatin-accessibility QTL and expression QTL
colocalize with largely distinct risk loci. It supplies cell-type regulatory
context that tissue-level eQTL analysis cannot provide (Elison et al., 2025).
These external studies extend the assay boundary that the Gene Catalog records.
They are credited source comparisons, not results of this paper.

The gene and module results of Kumak et al. are likewise source findings, not
results of this paper. Of their NAS-associated and fibrosis-associated genes,
944 and 1,013 map to the Catalog, and 258 and 301 of these reach the canonical
bulk call here. Of those present in both the bulk and colocalization tables, 661
of 902 and 714 of 943 had no multi-signal colocalization result. The Catalog
keeps this state separate from a tested negative. Among the 241 NAS and 228
fibrosis genes that were evaluated, 36 and 28 colocalized, against 28.9 and 24.9
expected for genes of the same expression level (Mantel-Haenszel OR 1.32, 95% CI
0.91 to 1.92, and OR 1.16, 95% CI 0.76 to 1.76; Holm-adjusted *p* 0.29 and 0.51,
Holm over the two lists). Thus, among genes evaluated by multi-signal COLOC, the
colocalization rate of these genes did not differ detectably from that of genes
of the same expression level (power in Supplementary Note 15). Kumak et al.
prioritize genes by combining meta-analysis evidence with centrality in a
curated protein interaction network. The 117 programs and 43 systems used here
were fixed before any outcome was read and are never reweighted by network
position. A sparsely annotated gene therefore keeps its assigned weight, and an
inadequately observed gene is recorded as untestable rather than ranked low.

The stage-stratified analyses add biological context without establishing
change over time. Adjacent-stage contrasts and cell-resolved programs describe
stage-ordered remodeling across different participants. The externally anchored
molecular continuum, adapted from Kamzolas et al. (2026), adds a separate
molecular view without redefining disease. It aligns strongly with the
fibrosis-based pathway, NMF and Hotspot results, which is consistent with the
major biological programs being preserved. This alignment was computed on the
same data and is not independent validation. Stage-adjusted associations show
that molecular variation within histologic bins tracks the continuum, and which
programs carry that variation depends on adjustment for cell-type composition.
Histology therefore remains the interpretable clinical anchor, and the continuum
adds molecular resolution within a stage. All 117 Hotspot program definitions
and 43 system memberships remain fixed. An outcome-blind map extends this
comparison across complete pathway and Hotspot layers. Stage and continuum
effects align broadly, but only a minority of systems keep one sliding-window
direction across all five cohorts. This supports shared but cohort-variable
remodeling, not a universal molecular sequence or ontology. Program-pair
covariance is reported in Results. It could not be tested in GSE213621, where
the program-group permutation control failed. Bulk covariance also mixes cell
abundance with program activity, so it does not show coordination between cell
types. Neither analysis identifies transitions within individuals, a dominant
fibrosis switch or therapeutic windows. Those claims need repeated sampling of
the same participants and perturbation experiments.

For users of the Resource, the main lesson of Figure 7 is that each model score
must be read beside its exposure audit and its task-native baseline.
Our per-cohort exposure audit showed that all 109 GSE268273
validation samples were in the BulkFormer pretraining data, so that model could not
be evaluated on samples it had not seen. Against task-native baselines, the paired
difference between the severity instrument and the in-cohort ridge has a 95%
interval that includes zero, and GSE268273 was also used to choose the
expression transform. A per-region model of H3K27ac that adds the expression of genes near each
region (the cis head) scored higher than a model on global expression components
alone, but its gain over a reduced-rank RNA model like that of Figure 7B was
smaller (Supplementary Note 14).

There are several limitations to this study. The liver-eQTL panel is
reused and mostly European, so some non-European GWAS lack ancestry-matched
molecular QTL data. The regulatory-genetic layer is based on expression QTL, so
protein-altering and protein-level mechanisms that leave the transcript
unchanged, such as MTARC1, fall outside its resolution. This layer is also
tissue-level, not cell-type-resolved. Source-reported protein-QTL and
single-cell chromatin studies supply distinct assay contexts, but no single
layer is a complete causal filter. Most genetically supported genes come from
liver-enzyme rather than direct MASLD traits. Partitioned polygenic scores
already show that MASLD risk variants do not form one class (Jamialahmadi et
al., 2024). The direct-versus-enzyme split used here is a different, trait-based
partition. The class supported on both axes is small, and every context assay
has source-specific coverage. No accepted spatial analysis provides a
MASLD-versus-healthy disease-effect estimand, and this paper contains no new
experimental perturbation. The GTEx comparison measures histologic fat of
unrecorded cause in postmortem liver. Its slide tissue-area negative control
missed its prediction (Spearman *rho* 0.354, *p* = 0.0564) and came close to its
stop rule. Ischemic time was also correlated with the same genes' disease
effects, with the opposite sign and less strongly than fat (*rho* −0.350 against
0.636). The caQTL reanalysis rests on one 138-donor map. Its chr1 and chr19
positive controls were seen in a pilot and are not blind, and at least 45 of the
191 colocalized signals were already reported by its authors. We did not check
overlap with published three-way eQTL, caQTL and GWAS colocalizations
(Broadaway et al., 2024). We record these as Catalog limitations and do not
infer mechanisms from source-owned data.

Further limitations concern units of analysis and models. Genotypes called from RNA-seq
reads identified repeated libraries and conflicting labels, but these calls,
sample-key checks and expression-based QC do not give complete genotype-based
identity confirmation, which was not established here. The bulk units therefore remain genetic-identity
clusters, not verified donors. Gene-level totals can hide opposing isoform
changes. Differential gene expression concerns total gene abundance,
differential transcript expression concerns individual transcript abundance, and
differential transcript usage concerns relative isoform proportions within a
gene. This Resource tests the first, not isoform switching. Short reads shared
among isoforms add assignment uncertainty, and transcript estimates depend on
the annotation. Oakes and Pachter emphasized these distinctions and the need to
separate demonstrated responses from exploratory physiological explanations
[ref: Oakes and Pachter, doi:10.1186/s12863-026-01440-8]. Model validation is specific to its task, cohort and
training exposure, and rank correlation is not absolute-stage calibration. The
GSE193066 repeat-biopsy intraclass correlation (ICC) and SD 0.6433 use
source-declared pairs whose identity was not qualified. They cannot establish
primary within-person reliability or label-error calibration, and ceiling values
derived from them remain sensitivity scenarios. Serial biopsies can include
shared technical effects and true change, and uncertainty conditional on fitted
predictions omits retraining variation. The comparator models were fitted
within each test cohort, so they had already seen that cohort's data.
Comparison with them therefore cannot show that the models perform equally on
unseen cohorts, or that one model class is better in general.

Two questions stay open. We do not assign genes or programs separately to steatosis, inflammation or fibrosis; fibrosis stage is the severity anchor, and NAS enters only as a sum or as an adjustment set. Colocalization here tests single variants, so the Resource does not test whether combinations of variants act differently, and the local haplotype at GNMT is unresolved. Both questions are recorded as indeterminate, not negative.

Within these limits, the Resource gives each gene an experiment matched to its evidence. Genetically anchored candidates motivate allele-to-gene and
context-appropriate perturbation without requiring disease-state expression.
Disease-state-associated candidates motivate tests of active tissue function and
reversibility without being labeled inherited drivers. Candidates supported on
both axes motivate experiments that connect the two relationships in one model.
Untestable or source-dependent candidates first need the missing decisive
measurement. The Gene Catalog and the source tables needed to check each call
are deposited with this paper (Data availability), so that others can look up
these distinctions for any gene instead of relying on a single aggregate rank.

# References

1. Kamzolas, I. et al. A data-driven framework reconstructs the molecular continuum of human MASLD progression. *Nat Metab* **8**, 1545–1562 (2026). https://doi.org/10.1038/s42255-026-01543-7

2. Rivas, M. A. et al. Effect of predicted protein-truncating genetic variants on the human transcriptome. *Science* **348**, 666–669 (2015). https://doi.org/10.1126/science.1261877

3. Kendall, T. J. et al. An integrated gene-to-outcome multimodal database for metabolic dysfunction-associated steatotic liver disease. *Nat Med* **29**, 2939–2953 (2023). https://doi.org/10.1038/s41591-023-02602-2

4. Alvarez-Sola, G. et al. LiRNA: An interactive atlas of human liver RNAseq databases. *JHEP Rep* **8**, 101914 (2026). https://doi.org/10.1016/j.jhepr.2026.101914

5. Piras, I. S. & DiStefano, J. K. Comprehensive meta-analysis reveals distinct gene expression signatures of MASLD progression. *Life Sci Alliance* **7**, e202302517 (2024). https://doi.org/10.26508/lsa.202302517

6. Gobeil, É. et al. Integrative genetic and liver transcriptomic analyses identify TRIB1AL as a target for steatotic liver disease. *J Clin Endocrinol Metab* **111**, 1980–1993 (2026). https://doi.org/10.1210/clinem/dgag043

7. Saliba-Gustafsson, P. et al. A functional genomic framework to elucidate novel causal metabolic dysfunction-associated fatty liver disease genes. *Hepatology* **82**, 165–183 (2025). https://doi.org/10.1097/HEP.0000000000001066

8. Jamialahmadi, O. et al. Partitioned polygenic risk scores identify distinct types of metabolic dysfunction-associated steatotic liver disease. *Nat Med* **30**, 3614–3623 (2024). https://doi.org/10.1038/s41591-024-03284-0

9. Currin, K. W. et al. Genetic effects on chromatin accessibility uncover mechanisms of liver gene regulation and quantitative traits. *Genome Res* **35**, 1485–1502 (2025). https://doi.org/10.1101/gr.279741.124

10. Broadaway, K. A. et al. Liver eQTL meta-analysis illuminates potential molecular mechanisms of cardiometabolic traits. *Am J Hum Genet* **111**, 1899–1913 (2024). https://doi.org/10.1016/j.ajhg.2024.07.017

11. Ghodsian, N. et al. Electronic health record-based genome-wide meta-analysis provides insights on the genetic architecture of non-alcoholic fatty liver disease. *Cell Rep Med* **2**, 100437 (2021). https://doi.org/10.1016/j.xcrm.2021.100437

12. Ji, Y. et al. Human liver methionine cycle: MAT1A and GNMT gene resequencing, functional genomics, and hepatic genotype-phenotype correlation. *Drug Metab Dispos* **40**, 1984–1992 (2012). https://doi.org/10.1124/dmd.112.046953

13. Gobeil, É. et al. An integrated proteogenomic investigation of the human liver uncovers molecular drivers of steatotic liver disease. Preprint at medRxiv https://doi.org/10.64898/2026.06.04.26354903 (2026).

14. Kumak, E., Darde, T. & Konu, O. Network-based meta-analysis maps stage-dependent molecular programs in MASLD through MASLD-META NETWORK application. Preprint at bioRxiv https://doi.org/10.64898/2026.08.26.747338 (2026).

15. Filliol, A. et al. Opposing roles of hepatic stellate cell subpopulations in hepatocarcinogenesis. *Nature* **610**, 356–365 (2022). https://doi.org/10.1038/s41586-022-05289-6

16. Guilliams, M. et al. Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches. *Cell* **185**, 379–396.e38 (2022). https://doi.org/10.1016/j.cell.2021.12.018

17. Vu, H. et al. Progressive fibrosis in human MASLD is associated with spatially linked transcriptomic signatures of metabolic reprogramming and senescence. *JHEP Rep* **8**, 101657 (2026). https://doi.org/10.1016/j.jhepr.2025.101657

18. Boesch, M. et al. Integrated multi-omics identifies distinct macrophage alterations during progression of metabolic dysfunction-associated steatohepatitis. *Nat Genet* **58**, 1295–1307 (2026). https://doi.org/10.1038/s41588-026-02600-3

19. Yakubovsky, O. et al. A spatial atlas of the healthy human liver from live donors. *Nature* **653**, 1148–1157 (2026). https://doi.org/10.1038/s41586-026-10377-y

20. Weiss, C. A. M. et al. Single-cell spatial proteomics maps human liver zonation patterns and their vulnerability to disruption in tissue architecture. *Nat Metab* **8**, 741–756 (2026). https://doi.org/10.1038/s42255-026-01459-2

21. Hoang, S. A. et al. Gene expression predicts histological severity and reveals distinct molecular profiles of nonalcoholic fatty liver disease. *Sci Rep* **9**, 12541 (2019). https://doi.org/10.1038/s41598-019-48746-5

22. Govaere, O. et al. Transcriptomic profiling across the nonalcoholic fatty liver disease spectrum reveals gene signatures for steatohepatitis and fibrosis. *Sci Transl Med* **12**, eaba4448 (2020). https://doi.org/10.1126/scitranslmed.aba4448

23. Pantano, L. et al. Molecular characterization and cell type composition deconvolution of fibrosis in NAFLD. *Sci Rep* **11**, 18045 (2021). https://doi.org/10.1038/s41598-021-96966-5

24. Ramachandran, P. et al. Resolving the fibrotic niche of human liver cirrhosis at single-cell level. *Nature* **575**, 512–518 (2019). https://doi.org/10.1038/s41586-019-1631-3

25. Elison, W. et al. Single cell multiomics reveals drivers of metabolic dysfunction-associated steatohepatitis. Preprint at medRxiv https://doi.org/10.1101/2025.05.09.25327043 (2025).

26. Sveinbjornsson, G. et al. Multiomics study of nonalcoholic fatty liver disease. *Nat Genet* **54**, 1652–1663 (2022). https://doi.org/10.1038/s41588-022-01199-5

27. Oakes, C. & Pachter, L. Transcriptomic responses to endurance exercise training in rats. *BMC Genom Data* **27**, 74 (2026). https://doi.org/10.1186/s12863-026-01440-8

28. Gribben, C. et al. Acquisition of epithelial plasticity in human chronic liver disease. *Nature* **630**, 166–173 (2024). https://doi.org/10.1038/s41586-024-07465-2

29. Andrews, T. S. et al. Single-cell, single-nucleus, and spatial RNA sequencing of the human liver identifies cholangiocyte and mesenchymal heterogeneity. *Hepatol Commun* **6**, 821–840 (2022). https://doi.org/10.1002/hep4.1854

30. Namjou, B. et al. GWAS and enrichment analyses of non-alcoholic fatty liver disease identify new trait-associated genes and pathways across eMERGE Network. *BMC Med* **17**, 135 (2019). https://doi.org/10.1186/s12916-019-1364-z

31. Anstee, Q. M. et al. Genome-wide association study of non-alcoholic fatty liver and steatohepatitis in a histologically characterised cohort. *J Hepatol* **73**, 505–515 (2020). https://doi.org/10.1016/j.jhep.2020.04.003

32. Liu, Y. et al. Genetic architecture of 11 organ traits derived from abdominal MRI using deep learning. *eLife* **10**, e65554 (2021). https://doi.org/10.7554/eLife.65554

33. Haas, M. E. et al. Machine learning enables new insights into genetic contributions to liver fat accumulation. *Cell Genom* **1**, 100066 (2021). https://doi.org/10.1016/j.xgen.2021.100066

34. van der Meer, D. et al. The link between liver fat and cardiometabolic diseases is highlighted by genome-wide association study of MRI-derived measures of body composition. *Commun Biol* **5**, 1271 (2022). https://doi.org/10.1038/s42003-022-04237-4

35. Sinnott-Armstrong, N. et al. Genetics of 35 blood and urine biomarkers in the UK Biobank. *Nat Genet* **53**, 185–194 (2021). https://doi.org/10.1038/s41588-020-00757-z

36. Kanai, M. et al. Genetic analysis of quantitative traits in the Japanese population links cell types to complex human diseases. *Nat Genet* **50**, 390–400 (2018). https://doi.org/10.1038/s41588-018-0047-6

37. McCarthy, D. J. & Smyth, G. K. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* **25**, 765–771 (2009). https://doi.org/10.1093/bioinformatics/btp053

38. Pinheiro, J., Bates, D. & R Core Team. nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-168. https://CRAN.R-project.org/package=nlme (doi:10.32614/CRAN.package.nlme).

39. Gerard, D. Data-based RNA-seq simulations by binomial thinning. *BMC Bioinformatics* **21**, 206 (2020). https://doi.org/10.1186/s12859-020-3450-9

40. Zhu, B. et al. Integrative analyses elucidate transcriptional regulatory functions of risk alleles for metabolic liver disease. *Nat Genet* **58**, 1280–1294 (2026). https://doi.org/10.1038/s41588-026-02617-8

41. Myojin, Y. et al. Genome-wide CRISPR screen identifies ACSL3 as a regulator of lipotoxicity and progression of MASLD. *Hepatol Commun* **10**, e0884 (2026). https://doi.org/10.1097/HC9.0000000000000884

# Figure legends

## Figure 1

**Figure 1. A multimodal human MASLD Resource keeps inherited liver-trait evidence separate from disease-state RNA.**

**(A)** Resource arms, each shown with its biological unit, phenotype provenance, source reuse and assay applicability. Genetics: 35 Tier-1/2 GWAS strata, direct disease traits kept apart from liver-enzyme proxies. Bulk RNA: five cohorts, 774 genetic-identity clusters (138 control, 636 disease), not verified donors. Single cell: 7 datasets, 104 source donors, 102 analyzed; the inferential n differs by contrast. Protein, chromatin and spatial data add context in native units. All arms feed the MASLD Gene Catalog (28,233 gene identities). Reuse of a public source is not independent validation.

**(B)** Overlap of genes colocalized with liver eQTL (SuSiE PP.H4 > 0.5) and canonical disease-versus-control DEGs (BH q < 0.05 across 24,115 genes, absolute log2FC > 0.5; 774 genetic-identity clusters), within the 3,769 genes with a multi-signal COLOC state. Counts are genes: 426 colocalized, 220 DEGs, 26 in both. Gray bars show the exact null of the overlap within 10 AveExpr bins; darker bars mark its central 95% (13 to 30). Dashed line, expected overlap (21.2); filled point, observed (26); open circle, expected without stratification (24.9); gray diamond, 13 DEGs among 426 expression-matched non-colocalized genes (exact p = 0.052), one control count rather than a second expectation. Lower row: Mantel-Haenszel OR 1.29 (95% CI 0.84 to 2.01; exact stratified two-sided p = 0.28; one prespecified test, no correction), log axis. Expression is not fully controlled; no direction is read. Outside the comparison: 935 canonical DEGs with only a single-signal ABF result and 699 with no tier-1/2 COLOC record (503 lncRNA); these are untested, not negatives.

**(C)** Trait-class support of the 436 genes colocalized with liver eQTL (SuSiE PP.H4 > 0.5 in any of 35 main tier-1/2 GWAS). Each bar is one trait class: ALT/AST/GGT (20 studies, 424 genes) and NAFLD/NASH/PDFF (15 studies, 22 genes); ten genes are in both. Segments give each gene's SuSiE state in the other class: also supported (10 each), tested at PP.H4 ≤ 0.5 (46 enzyme, 11 direct) or not SuSiE-tested (outlined; 368 enzyme, 1 direct). Not SuSiE-tested means untested, not negative. Testing was unequal: 3,888 genes tested in at least one enzyme study, 392 in at least one direct study. Among 362 genes tested in both, 56 are enzyme-supported and 21 direct-supported. Counts are genes.

**(D)** Analysis flow across Figures 1 to 7. The same source-qualification rules (phenotype, biological unit, coverage and applicability, source dependence) apply to every analysis, and each reports a machine-readable evidence state, the unresolved alternative and the next discriminating experiment. Inherited liver-trait evidence (Figure 2), multicohort disease-state remodeling (Figure 3), donor-level programs (Figure 4) and molecular and tissue context (Figure 5) lead to the Gene Catalog (Figure 6) and model evaluation (Figure 7). Arrows show analysis order, not causal links.

## Figure 2

**Figure 2. Liver-trait GWAS colocalize with liver eQTL at 436 genes, mostly through liver-enzyme traits.**

**(A)** Three unlinked counts over 35 Tier-1/2 GWAS strata: studies by ancestry; fine-mapped loci (360 trait-specific, 265 unique) by ancestry of top PIP, with 54 tied loci as a separate segment; and 996 colocalized genes (436 multi-signal, 560 single-signal only), each counted once in the band of its largest PP.H4.

**(B)** Largest multi-signal PP.H4 per gene against SuSiE PIP of the same lead variant in the same GWAS (n = 866 genes; 327 colored, 539 gray outside the display filter). Shape gives lead consequence. Genes whose largest posterior comes from the unfine-mapped MVP European strata have no PIP and are not drawn; absence is not a negative call.

**(C)** Consequence class of one lead COLOC variant per gene (996 genes, 933 variants), first matching hg19 knownGene class in the order coding or splice, UTR, promoter, intron, intergenic; coding includes synonymous. A gene-level summary, not the posterior composition of credible sets.

**(D)** Fine-mapping posterior in four exclusive consequence classes, averaged over 1-Mb clusters for direct MASLD/liver-fat traits (207 credible sets, 19 clusters) and liver-enzyme traits (6,208 sets, 228 clusters). Inset: direct minus enzyme protein-altering plus splice probability, 19.5 percentage points (95% cluster bootstrap 6.9 to 34.0). Ancestry differs between arms, so this is not a trait-class contrast.

**(D, ancestry rows)** The same difference by GWAS ancestry. Thick black bars, 95%; thin gray bars, 98.3% cluster-bootstrap intervals (Bonferroni over EUR+AMR, EUR and AMR). Pooled and all-non-EUR rows are descriptive. No member receives a trait-class call.

**(E)** Ancestry of the GWAS supporting each of the 436 multi-signal genes (SuSiE PP.H4 > 0.5); a per-gene, not per-locus, partition. Of 340 EUR-only genes, 127 were tested without support in non-EUR studies and 213 were untested, which is not the same as unsupported. Black strip, non-EUR ABF PP.H4 > 0.5. The eQTL panel is European, so non-EUR GWAS are not ancestry-matched.

**(F)** All 13 genes with PP.H4 > 0.5 (SuSiE or ABF) in four of five ancestry groups. Circle area, maximum PP.H4 over each group's studies. Dark blue, SuSiE; light blue, ABF only; gray, neither above 0.5. Asterisk, four groups by SuSiE alone. These partitions do not establish biological differences between ancestry groups.

**(G)** RORA in MVP NAFLD EUR (N = 432,709; 29,342 cases). Lead p = 7.6e-8; no fine-mapping record, so no PIP track. SuSiE PP.H4 = 0.983.

**(H)** FABP1 in UKBB ALT (N = 343,850). Lead 2:88423510 (PIP 0.57) and COLOC top variant 2:88424066 (PIP 0.43) are hollow diamonds. SuSiE PP.H4 = 0.998. FABP1 is not differentially expressed in bulk liver (log2FC −0.04, padj 0.51).

**(I)** Key for (G) and (H). Tracks: GWAS p colored by r² to the lead (EUR LD), SuSiE PIP, Broadaway liver eQTL (N = 1,183) [ref: Broadaway et al., liver eQTL] and genes (500 kb, hg19). Gene fill is bulk log2FC, MASLD versus control, over 774 genetic-identity clusters (not verified donors), capped at ±0.3; gray, padj ≥ 0.05; outlined white, outside the bulk gene universe.

See also Figure S2.

## Figure S2

**Figure S2. Fine-mapping and colocalization detail, related to Figure 2.**

**(A)** The 35 prespecified Tier-1/2 GWAS strata, in 33 rows colored by GWAS ancestry (Pan-UKBB Central/South Asian labeled SAS). Left, loci with a SuSiE or CARMA credible set (log scale); the MVP European ALT, AST and NAFLD strata were not fine-mapped. Right, genes with multi-signal (SuSiE-coloc) or single-signal (coloc.abf) PP.H4 > 0.5. A gene colocalized in several strata counts in each row; the union over the 35 strata is 996 genes.

**(B)** Independent fine-mapped loci. Credible-set windows (1,665 in the 32 fine-mapped strata) are padded by ±250 kb and merged across strata, giving 265 loci (182 in one stratum, 83 in two or more). Method bars are split by best variant PIP. SuSiE and CARMA are counted on the merged loci; SuSiEx and meSuSiE on their own MVP cross-ancestry loci, a different locus set. Without padding there are 289 loci.

**(C)** The 11 merged loci of (B) with a result from all four methods (cross-ancestry loci matched within 500 kb), counted by best PIP. SuSiE and CARMA values are maxima over single-ancestry strata; SuSiEx and meSuSiE values come from one joint MVP fit, so the comparison is not like-for-like. One locus is 9 percentage points, and the panel does not support a difference between methods.

**(D)** EFHD1 locus for AST on chromosome 2 (hg19). Tracks show European (UKBB) and East Asian (BBJ) GWAS association colored by r² to the shared lead, BBJ SuSiE PIP (filled, credible-set member; UKBB AST was not fine-mapped here) and the EFHD1 liver eQTL. Genes are filled by bulk log2FC when they have a canonical call (EFHD1 −1.10) and gray otherwise.

**(F)** Method agreement in the 35 strata. Left, SuSiE-coloc against coloc.abf PP.H4 for the 12,892 gene-stratum pairs with both values (Pearson r = 0.74). Printed coloc.abf counts are over all 664,022 coloc.abf pairs. Right, SuSiE against CARMA PIP for 293,700 variants in the 13 strata where CARMA was run (r = 0.76).

**(G)** Genes colocalized at PP.H4 > 0.5, > 0.8 and > 0.9 by each method, split by ancestry support: European if any European stratum passes, non-European-only otherwise. The 32 SuSiE-coloc non-European-only genes at PP.H4 > 0.5 match Figure 2E.

**(H)** The 244 genes with SuSiE-coloc PP.H4 > 0.9, by chromosome (none on chromosome 18). Columns give the maximum PP.H4 in each GWAS ancestry; gray cells have no COLOC result. Side strips give the class of the lead variant (top SNP of the gene's best multi-signal result) and the bulk call; white marks genes outside the bulk universe (not tested). The 32 labeled genes pass 0.9 in two or more ancestries.

**(I)** The 436 SuSiE-coloc genes (PP.H4 > 0.5) by genome position and best PP.H4 (dashed lines at 0.8 and 0.9). Color gives lead-variant class as in (H); point size gives the number of ancestries with PP.H4 > 0.5. Labeled genes pass 0.9 in two or more ancestries or have a coding or splice lead; 144 other genes at PP.H4 ≥ 0.95 are not labeled.

**(J)** Gene maximum PIP from one joint MVP meSuSiE fit for ALT, AST and NAFLD, at 13 known liver-trait loci (TM6SF2 and GATAD2A share a row). Open circles, European-specific component; filled circles, maximum over the 11 shared components. Neither is a European-only fine-mapping. Genes are linked to variants through liver cis-eQTL pairs (Broadaway et al.), so a maximum can lie outside the gene body. Nine of the 79 loci, including PNPLA3, TRIB1, RECQL4 and SERPINA1, have a European-component maximum of 0.4217 to 0.4218, which likely reflects how meSuSiE divides posterior mass, not a separate European estimate.

**(K)** Annotation class of the 1,263 unique COLOC top variants from the 1,753 gene-stratum pairs with PP.H4 > 0.5 by either method (996 genes), each variant counted once and assigned by priority (coding first, intergenic last; promoter −2 kb to +200 bp). The unit is the variant, not the one lead per gene of Figure 2C.

**(L)** F2RL1 PP.H4 by method in European (UKBB) and South Asian (Pan-UKBB) GGT GWAS, beside bulk log2FC (+0.78, 95% CI 0.62 to 0.94; BH padj = 1.1 x 10^-18). South Asian GGT has no multi-signal result, so that cell is untestable, not negative.

**(M)** The 256 meSuSiE credible sets in the 79 MVP ALT, AST and NAFLD loci, by ancestry sharing. East Asian GWAS are in the fit, but no set is specific to East Asian ancestry, so there is no East Asian bar.

**(O)** Six genes compared between liver eQTL colocalization (this study; the larger of SuSiE-coloc and coloc.abf PP.H4 over the 35 strata) and protein-QTL colocalization (Gobeil et al., 2026b [ref: medRxiv doi:10.64898/2026.06.04.26354903]). Row text gives the eQTL support state. HSD17B13 is below the 0.80 protein-QTL threshold used by Gobeil et al.; for GCKR the source reports opposite eQTL and pQTL effect signs.

**(P)** The six genes of (O) as eQTL PP.H4 against protein-QTL PP.H4. Color gives eQTL support (multi-signal, single-signal only, or none above 0.5); the triangle marks the source-reported sign discordance (GCKR). Dashed lines mark eQTL PP.H4 = 0.5 and protein-QTL PP.H4 = 0.8.

**(Q)** 17 genes from a published liver single-cell chromatin study (Elison et al. [ref: Elison et al., liver single-cell chromatin and caQTL study]), grouped by the source's layer. Columns give the maximum eQTL PP.H4 (either method) over the 15 direct MASLD/MASH/liver-fat strata and over the 20 liver-enzyme strata, the published finding and cell type, and the bulk call with BH padj. An asterisk marks PP.H4 > 0.5 from coloc.abf only (multi-signal PP.H4 at most 0.5 or not run).

**(R)** 27 genes whose locus carries an MPRA differential-activity variant in a published liver MPRA (Zhu et al. [ref: Zhu et al., liver MPRA of disease-associated variants]). Columns as in (Q), plus the strongest allelic MPRA log2FC in HepG2 or LX-2 cells and every map in this study with support. NCAN and DLG5 are outside the 24,115-gene bulk universe and are not evaluable, which differs from no call.

**(S)** Coding plus canonical-splice share of fine-mapped posterior (Figure 2D), split within the direct class. PDFF (three European GWAS) and NAFLD plus NASH are compared with all enzyme traits and with each other; PDFF is also compared with European enzyme sets only, because all PDFF sets are European. Bars are 95% physical-cluster bootstrap intervals (10,000 draws; 20,000 with subclasses resampled independently for PDFF minus NAFLD plus NASH).

**(T)** Phenotype provenance of the 35 strata: liver enzyme (ALT, AST or GGT; 20 strata, 5 European), direct MASLD/MASH diagnosis (12; 9 European) and MRI-PDFF or histologic steatosis (3, all European). Non-European GWAS are paired with the European eQTL panel, which limits colocalization.

**(U)** SuSiE-coloc genes (PP.H4 > 0.5) per trait-ancestry cell, and those with a lead within ±250 kb of a fine-mapped window from the same study; the Windows column counts fine-mapped windows. Totals are 436 genes, 100 with such a lead, and 1,665 windows. European enzyme cells have few such genes because MVP European strata were not fine-mapped and UKBB strata only at few loci.

**(V)** Regulatory context of noncoding fine-mapped posterior after re-annotation by position (Supplementary Note 6), as a share of noncoding posterior summed over all credible sets of a class. Contexts overlap and are not summed. Pooled values mix ancestries. Within AMR, only the no-context difference passed its calibration control, and no difference was detected.

**(W)** Target-gene biotype of the 436 primary SuSiE genes: 424 protein-coding and 12 lncRNA. Targets are gene-level and do not resolve transcript isoforms.

Bulk calls use the 774-unit disease-versus-control fit (genetic-identity clusters, not verified donors). A canonical call is BH padj < 0.05 over 24,115 genes and |log2FC| > 0.5.

## Figure 3

**Figure 3. Bulk liver RNA from five cohorts defines disease- and fibrosis stage-associated genes.**

Units are 774 genetic-identity clusters (138 control, 636 disease), not verified donors. The pooled fit tests 24,115 genes and gives 1,854 canonical calls (BH padj < 0.05 and |log2FC| > 0.5).

**(A)** Metadata recorded by each cohort (tiles) and units (bars, disease and control). GSE162694 records NAS for 117 and fibrosis for 109 of 140 units; the 31 without fibrosis are controls.

**(B)** Units by exact fibrosis stage and NAS (364, controls included). The top row counts all 387 units with an exact stage (F0/F1/F2/F3/F4 = 103/102/88/69/25), used in (E) and (F). GSE213621 records only fibrosis bins.

**(C)** Per-cohort calls (BH FDR < 0.05 and |log2FC| > 0.5 within each cohort's 24,115-gene family). Top, calls by direction. Middle, genes per overlap class in the 13,761-gene union; darker bars include a GSE126848 call. Bottom, mixed-sign share of multi-cohort genes, with and without GSE126848. Shared membership is not replication. GSE126848 calls (7 disease, 23 control clusters) also track library depth, and GSE135251 calls cannot be matched for genetic ancestry (6 of 10 controls of African genetic similarity).

**(D)** Group-blind PCA of 113 strict controls and 279 advanced-disease clusters (NAS ≥ 5 or fibrosis ≥ F3) in four cohorts, adjusted for cohort and inferred sex. Open gray circles, controls; filled circles, exact fibrosis stage; open triangles, 77 GSE213621 units in the F3–F4 bin. Top, PC1 density by entry rule. These units overlap the pooled fit, so this is a sensitivity analysis, not independent replication.

**(E)** Top, genes differing from F0 at each stage (BH q < 0.05 within each 24,115-gene family, no fold-change filter; 1,338 / 4,266 / 6,190 / 7,043 for F1–F4), with cohort and inferred-sex terms. Lighter inner bars, genes also with |log2FC| > 0.5. F0 holds 14 source controls and 89 disease biopsies. The F1 count rests mainly on GSE130970; without it, 14 genes remain. Bottom, F0-centered activity of five k6 NMF axes (continuous, not patient classes) with 95% intervals (estimate ± 1.96 SE). F0-reference and cohort-omission views are in Figure S3D–F.

**(F)** Fifty-four genes of the 436-gene primary SuSiE set with BH FDR < 0.05 in at least one adjacent fibrosis contrast (205 / 190 / 157 / 94 clusters). Color, higher-minus-lower stage log2FC, clipped at ±0.688 (triangles). Filled dots, BH q < 0.0125 (5% split across four contrasts; 15 genes); open dots, 0.0125 ≤ q < 0.05. FDR is controlled within each contrast, not over their union. No F0→F1 mark keeps BH q < 0.05 without GSE130970, and no F1→F2 mark without GSE162694. Right, PP.H4 and GWAS trait of each gene's strongest tier-1/2 driver. Selection is conditioned on stage significance; no enrichment or causal claim is made.

**(G)** Disease-state |log2FC| against unsigned SuSiE COLOC PP.H4 for 3,769 genes with a SuSiE test. Dotted lines mark the class rules. The lower strip expands PP.H4 < 0.05, where 2,939 Neither genes are binned. Open circles, seven SuSiE-untestable genes, not counted as Neither. Counts use the native RNA gene universe, not the Gene Catalog universe (Figure 6). Fisher and Spearman tests are nominal and unstratified; the prespecified stratified test is in Figure 1B. Inset, lncRNA versus protein-coding: odds of a call over 21,501 genes, unadjusted (open) and adjusted for expression and detection (filled), 95% CI; share of calls up; conditional odds of up in 520 covariate-matched pairs (Wald 95% CI).

See also Figure S3.

## Figure S3

**Figure S3. Bulk transcriptomic sensitivities and continuum analyses, related to Figure 3.**

Except in (H) and (X), units are drawn from the 774 genetic-identity clusters, which are not verified independent donors. GSE213621 records only Control, F0–F1, F2 and F3–F4, so stage displays show it apart from cohorts with exact stages. Category means and overlapping windows are descriptive.

**(A)** Principal component analysis of the 774 clusters (138 control, 636 disease; five cohorts). Input is voom log-CPM of the 2,000 most variable genes after regressing out cohort and inferred sex, without protecting the disease label. PC1 explains 7.9% and PC2 5.2% of variance. Gray open circles are controls and magenta points are disease. The top strip shows PC1 density per group.

**(B)** Log2FC for advanced disease (n = 279) versus strict control (n = 113) in four cohorts, against log2FC from the five-cohort disease-versus-control fit (774 units). Points are genes called in either contrast; dashed line, y = x. Over all 24,115 genes, Spearman rho = 0.953. Of 1,854 canonical calls, 1,755 are also called in the endpoint contrast (recovered, magenta) and 99 are not (subthreshold, pink). A further 1,195 are called only in the endpoint contrast (gray). No gene is called in opposite directions.

**(C)** Direction agreement of the 1,854 pooled calls in each cohort ranges from 78.3% to 100% (GSE126848, 84.36%). Leave-one-cohort-out fits on the fixed 24,115 genes give effect Spearman correlations of 0.870–0.995 and call retention of 70.7%–96.6%. Full refits that repeat filtering after cohort removal retain 66.7%–96.8% of calls with the same sign.

**(D)** Top, genes at BH q < 0.05 (no fold-change filter) for each fibrosis stage versus F0 (left) and each NAS group versus NAS0 (right), split into higher (above the axis) and lower (below) expression in the comparison group. Middle, F0- or NAS0-centered activity of the five available k6 NMF axes (P1–P5). Bottom, estimated hepatocyte and non-parenchymal fractions. The NAS0 reference (46 units) holds 45 source controls and 1 disease biopsy, so NAS 1–2, 3–4 and 5–8 versus NAS0 (5,338, 9,337 and 9,841 genes) mainly contrast disease with controls. Between grades, NAS 1–2 to 3–4 gives 3 genes (65 and 117 units) and NAS 3–4 to 5–8 gives 4,305 (117 and 167 units).

**(E)** F1, F2, F3 and F4 versus F0 gene counts (BH q < 0.05 in each 24,115-gene family) under three F0 references: all 103 F0 units (1,338 / 4,266 / 6,190 / 7,043), the 89 disease-only F0 units (783 / 3,946 / 6,459 / 6,929), and all 103 with a control-status term (814 / 3,679 / 6,166 / 6,640).

**(F)** F1 versus F0 gene counts (BH q < 0.05). Three-cohort model, 1,338 genes (205 units). Without GSE130970, GSE135251 or GSE162694, 14, 929 and 651 genes (155, 115 and 140 units). Each cohort alone, 252, 3 and 0 genes (50, 90 and 65 units).

**(G)** 426 one-to-one RNA-abundance-matched pairs of a genetically supported gene and a control gene. Interval-null calls (lfc = 0.25, BH q < 0.05) occur in the genetic member only in 37 pairs, the control only in 28, neither in 360 and both in one (paired OR = 1.321; nominal exact McNemar p = 0.321). This is a detectability sensitivity, not a test of risk-allele direction.

**(H)** Sample-level Spearman clustering in an earlier nine-cohort bulk sample set, after batch and sex correction that conditioned on disease status. Qualitative only: it does not use the five-cohort, 774-unit data and supplies no number reported in this study.

**(I–L)** The two fixed co-primary continuum scores against recorded fibrosis (I), in stage-adjusted anchor models (J), by stage distribution (K) and against each other (L; Pearson r = 0.999 in four cohorts and 0.979 in GSE126848, 30 units). Within-cohort percentile ranks are display coordinates, not disease probabilities. Stage displays use the four cohorts with recorded fibrosis. Inference is restricted to GSE162694 and GSE213621; source-overlap cohorts are descriptive. Stage-complete anchor models use 435 clusters, including the 68 GSE213621 controls, which have no recorded histology and are placed at stage 0.

**(M–R)** Continuum-associated transcripts: complete transcript effects (M), adjacent-stage comparisons (N), the fixed five-gene roster (O, P), all-gene disease and continuum concordance (Q), and projection of a disease signature derived from the source-overlap cohorts (R). Models adjust for fibrosis stage and inferred sex in the two evaluation cohorts. CYP2C19 and IGFBP7 use target-specific leave-one-gene-out axes. Window profiles in (P) show all five cohorts descriptively. These transcripts are not a redefined DEG family.

**(S–U)** All 50 Hallmarks (S), and six prespecified Hallmarks by continuum window (T) and fibrosis category (U), after removing the 143 observed genes of the 145-gene signature. BH covers the full declared family. Windows are nine 20-percentile-wide windows with 50% overlap, so each unit can appear in two windows; window displays use all five cohorts and fibrosis displays the four with recorded stage. Means and pointwise intervals do not establish support or order in time. EMT keeps 75.5% of its membership, below the 80% rule, and is labeled untestable with no curve drawn; the other five are testable.

**(V)** Fibrosis-stage trend β against stage-adjusted fixed-projection β for every testable pathway in Hallmark, KEGG, Reactome and the three GO collections, signature genes removed. Points are pathways, colored by whether they meet the stage criterion, the continuum criterion, both or neither. Units are the 435 stage-complete clusters of GSE162694 and GSE213621. Both axes use the same expression data, so this shows alignment, not independent validation or order in time within a person.

**(W)** −log10 complete-family q against adjusted meta-analysis β for every pathway in the six collections, for the fixed-projection score (left) and the cohort PC1 score (right).

**(X)** Residual rank (partial Spearman) correlations for each GNMT, MAT1A and CYP2C19 protein pair in the same 46 MASLD participants. Adjustments: baseline (acquisition batch, age, BMI and sex); baseline plus fibrosis stage; baseline plus NAS; baseline plus fibrosis stage and steatosis, ballooning and inflammation grades; baseline plus metabolic background (median z of three Hallmark metabolic sets). Bars are pointwise participant-bootstrap 95% intervals (2,000 resamples). The example is selection-conditioned, and adjustment does not establish mediation.

## Figure 4

**Figure 4. Single-cell programs show no corrected stage association but vary along a bulk molecular continuum.**

**(A)** The atlas integrates 275 libraries and 1,232,318 cells from 7 datasets (104 donors). We analyzed 273 libraries (1,232,285 cells) from 102 donors and defined 117 fixed Hotspot programs in 5 lineages.

**(B)** UMAP of 50,000 cells drawn uniformly from the 1,232,285 analyzed cells (102 donors), colored by cell type (analyzed cells per type in key). Dendritic cells include plasmacytoid and conventional subsets.

**(C)** Ordinal stage β for the 117 programs, ranked within lineage, from generalized least squares with per-dataset residual variance in 60 ambient-correctable donors (GSE202379, GSE244832; donors per lineage above each block). 17 of 53 T-cell donors have fewer than 20 T cells. Gray bars, β ± standard error times the within-dataset permutation critical |t| at 0.05/117; black bars, nominal 95% confidence intervals. No program passes the 117-program permutation BH and stability rule (smallest q = 0.323).

**(D)** Discovery against content-associated lineage for the 117 programs (donor-level same-atlas scores). Eighty match (shaded diagonal) and 37 do not. Labels mark the 7 calls that keep the same lineage in fewer than 95% of 1,000 donor bootstrap draws. Content association does not establish transcript origin or cell autonomy.

**(E)** Fibrosis effect (cohort-standardized program-score SD per recorded stage-code unit) for 112 eligible programs, ranked, in 364 genetic-identity clusters; NAS effects in Figure S4W. Filled points, random-effects meta-estimates colored by discovery lineage; gray bars, Knapp-Hartung 95% confidence intervals (t, 2 df); open marks, per-cohort estimates (circle GSE130970, n = 71; triangle GSE135251, n = 207; square GSE162694, n = 86). No association meets BH q < 0.05 over 224 tests (112 programs x fibrosis and NAS; smallest q 0.079 and 0.159), with or without composition adjustment. With z-based REML inference instead, 37 fibrosis and 14 NAS associations reach q < 0.05.

**(F)** Cohort-standardized Stromal ECM and Ductular injury scores in nine overlapping fixed-projection windows (774 genetic-identity clusters, five cohorts). Lines join window means; bands are ± 1 standard error; clusters per window in key. Dashed gray lines, cohorts overlapping continuum construction; solid lines, the two independent evaluation cohorts. The Stromal ECM score excludes 9 continuum-signature genes carrying 19.5% of its weight. Descriptive; window means do not imply order in time. Adjusted estimates in Figure S4V.

**(G)** Stage-adjusted continuum coefficient for each of 5,435 testable nodes (gray; pathways and Hotspot programs) in 43 fixed, signature-excluded systems, estimated in 367 stage-complete disease units (109 GSE162694, 258 GSE213621). Rows, systems by size (node count in brackets). Open circles, system medians with 95% bootstrap intervals (disease units resampled within cohort, B = 1,000). Four untestable nodes are counted at right.

**(H)** Fibrosis-stage coefficient against stage-adjusted continuum coefficient for the 43 systems, both donor-first in the same 367 units. Bars, 95% bootstrap intervals (B = 1,000). Color marks systems whose two intervals both exclude zero (19 both positive, 2 both negative, 1 opposite signs; 21 not resolved, gray). Spearman rho = 0.67 (bootstrap 95% interval 0.38 to 0.84); 34/43 share direction.

See also Figure S4.

## Figure S4

**Figure S4. Single-cell program, communication and continuum detail, related to Figure 4.**

Bulk panels use selected genetic-identity clusters, not verified independent donors. In continuum displays, dashed gray lines (†) mark GSE126848, GSE130970 and GSE135251, which overlap continuum construction and are descriptive only; solid lines mark GSE162694 and GSE213621. Windows are nine fixed-projection windows, 20 percentile points wide and centered at the 10th to 90th percentiles, so adjacent windows overlap. Window and category means are descriptive. Untestable is not a negative result.

**(A)** Hepatocyte native-score effect per stage step for Stromal ECM and Ductular injury in the 62-donor equal-variance sensitivity. Rows are the pooled fit, the two multi-stage datasets and five dataset omissions, with donor n on each row. GSE174748 and GSE185477 (Healthy only) and GSE189600 (one donor per stage) have no dataset row. Bars are nominal HC3 95% intervals, not new complete-family calls.

**(B)** MASH minus steatosis projected score in fibroblasts (Stromal ECM) and cholangiocytes (Ductular injury), raw and ambient-corrected, on the same cells and donors at the 50-cell threshold (GSE244832, n = 13 and 12 donors; GSE202379, n = 13 and 14). Bars are cohort-specific HC3 95% CIs. The two corrected tests form one BH family. Projected scores are not native Hotspot units.

**(C)** Median dataset-centered donor score of all 117 Hotspot programs in Healthy, Steatosis and Steatohepatitis donors, grouped by lineage. Right, GLS stage β from Figure 4C (60 ambient-correctable donors) with nominal 95% CI (black) and 117-program permutation interval (gray).

**(D)** Donor-collapsed communication scores for the fixed 13-pair roster across Healthy, Steatosis and Steatohepatitis (69 donors). Six pairs are testable, all with corrected q = 1. Seven are untestable because stage counts are insufficient; these are drawn at lower opacity and are not tested nulls.

**(E)** Adjusted associations of the two focal programs with the 13-pair roster (26-test BH family, twenty-two tests estimable). Two Ductular injury associations pass model-based and HC3 BH, each in 58 donors: hepatocyte→endothelial *CDH1*→*PTPRM* (beta = 0.4957, q = 3.628 x 10^-6, HC3 q = 1.167 x 10^-10) and endothelial→hepatocyte *ST6GAL1*→*EGFR*. Right, partial residuals for *CDH1*→*PTPRM*. These are RNA associations, not evidence of ligand secretion or receptor activation.

**(F)** Cohort-standardized Ductular injury and Stromal ECM scores across the nine windows in all five cohorts. Point size gives units per window.

**(G)** Mean Ductular injury and Stromal ECM scores (cohort z) by recorded fibrosis stage in GSE162694 (exact F0–F4; n = 35, 30, 27, 8 and 9) and by source category in GSE213621 (Control, F0–F1, F2 and F3–F4; n = 68, 87, 94 and 77). Bars are pointwise 95% intervals (mean ± 1.96 SD/√n). GSE162694 controls without a recorded stage (n = 31) are not plotted.

**(H)** All 117 Hotspot programs across the nine windows. Rows give cell type, module number and program name. Right columns give the adjusted meta-analysis effect and support status from the two evaluation cohorts.

**(I)** The same display restricted to the 93 of 117 programs (116 estimable) that are continuum-associated in the 435 stage-complete units of the two evaluation cohorts, including the 68 GSE213621 controls at stage 0.

**(J)** Fibrosis-stage trend β against stage-adjusted fixed-projection β for all 117 programs, colored by support: both (54), continuum only (39), stage only (10), neither (13); one program is not estimable. Panels J and L use the 435 stage-complete units, with the 68 GSE213621 controls placed at stage 0.

**(K)** Adjusted meta-analysis β of the ten fixed NMF axes (k4 P1–P4, k6 P1–P6) on the cohort PC1 and fixed-projection continuum scores, with 95% CIs. Per cohort (GSE162694, n = 109; GSE213621, n = 326), OLS of NMF score on continuum score, categorical fibrosis stage and sex; cohorts combined by fixed-effect inverse-variance meta-analysis. The 68 GSE213621 controls are included at stage 0.

**(L)** As in (J) for the ten NMF axes. k4-P2 and k6-P6 fail the fixed reconstruction rule for stage and have no stage-comparison point.

**(M)** Mean scores of the ten NMF axes (cohort z) in the nine windows in all five cohorts, including controls. Point size gives genetic-identity clusters per window (6 in GSE126848, 14 or 15 in GSE130970, 41 or 42 in GSE135251, 28 in GSE162694, 65 or 66 in GSE213621).

**(N)** Mean NMF axis scores (cohort z) by recorded fibrosis stage in the four cohorts that record it: GSE130970 (n = 23, 27, 9, 10 and 2 at F0–F4) and GSE135251 (n = 45, 45, 52, 51 and 14), dashed gray; GSE162694 (n = 35, 30, 27, 8 and 9), magenta; GSE213621 (Control, F0–F1, F2 and F3–F4; n = 68, 87, 94 and 77) on a separate axis, Control in gray. Bars are pointwise 95% intervals, omitted for n < 3.

**(O)** Signed regulon scores (cohort z) for HNF4A, RELA, SMAD3 and TCF7L2 across continuum windows in five cohorts. RELA and SMAD3 keep 74.3% and 75.8% of their original L1 weight and are untestable. Regulon scores are inferred, not direct TF activity.

**(P)** Paired-change β for Ductular injury in 54 source-declared GSE193066 biopsy pairs (108 libraries), outside the 774-unit data: change in program score regressed on change in continuum score and change in fibrosis stage. β = 0.644 (95% CI 0.447 to 0.840; q = 7.7 x 10^-7) on the fixed projection and 0.634 (0.442 to 0.826; q = 6.3 x 10^-7) on cohort PC1; BH over 117 programs. Pairs are not identity-verified, so this is not within-person support. Stromal ECM is untestable (73.8% of weight retained).

**(Q)** UMAP of gene-set membership (axes arbitrary) for 5,445 mapped features after signature-gene exclusion, colored by collection. Labels mark the 15 largest of 43 systems. Distances are not inferential.

**(R)** Fibrosis-stage β (score SD per stage) of each feature on the coordinates of (Q).

**(S)** Collection-balanced fibrosis-stage β against stage-adjusted continuum β for the 43 systems, colored by direction and sized by mapped features. In the 435 stage-complete units, including the 68 GSE213621 controls at stage 0, Spearman rho = 0.73 and 37 of 43 share sign; in the 367 disease units of Figure 4H, rho = 0.67 and 34 of 43 (donor-first). Descriptive, because gene sets are nested.

**(T)** Within-cohort collection-balanced system scores (z) across the nine windows in each of the five cohorts. Rows follow the order of (U); + and − give the sign of the late-minus-early contrast.

**(U)** Equal-cohort pooled heatmap of all 43 systems across the nine windows, a descriptive summary, not participant pooling. Rows are ordered by the pooled late-minus-early contrast (mean of windows at 80 and 90 minus mean at 10 and 20), for presentation only. Marks at right identify the 18 systems with strict maxT support (p < 0.05 in both cohorts on both co-primary axes, 435 units). These 18 come from models without composition terms; 14 lose the call in at least one composition-adjusted fit.

**(V)** Stage- and sex-adjusted and composition-adjusted continuum effects of Stromal ECM and Ductular injury in the 367 disease units (GSE162694, n = 109; GSE213621, n = 258), with composition adjustment under both cell-type fraction sets. Bars are HC3 95% CIs. No composition-adjusted program or system is named as a finding.

**(W)** NAS-axis companion to Figure 4E: program-score SD per recorded NAS point for 112 programs, from the same per-cohort OLS fits (sex, fibrosis stage, NAS; HC3) in 364 genetic-identity clusters, combined by REML random-effects meta-analysis with Knapp-Hartung 95% CIs. Colors give discovery lineage. The five largest positive effects are labeled. BH over the 224 tests of Figure 4E; no NAS association meets q < 0.05 (smallest q 0.159).

## Figure 5

**Figure 5. Protein, chromatin and spatial data qualify the genetic and RNA evidence in their native units.**

**(A)** Inputs and readouts. Each input set is fixed before its readout; no program is reselected.

**(B)** Liver protein against bulk mRNA change for 3,054 measured genes of a 10,044-gene prioritized universe. The universe is partly selected on RNA: rho is 0.487 in RNA-selected genes and 0.305 in the 664 genes selected on genetic evidence alone. Protein model: 46 MASLD versus 12 control participants, adjusted for acquisition batch, age, BMI and sex. Blue, protein-significant (BH p < 0.05 over 6,071 proteins) with the same sign as mRNA; magenta, opposite sign; gray, not significant. Bars split the 211 protein-significant genes by RNA call.

**(C)** 25 proteins fixed before the bulk refit, grouped by process. The same cohort helped select them, so the view is descriptive. Heatmap, protein z-scores, ordered by fibrosis stage then NAS. Middle, partial Spearman rho with each histology score within 46 MASLD participants, participant-bootstrap 95% intervals; filled, BH q < 0.05 over 125 tests. Smaller n are printed beside rows. Right, filled markers are calls (BH padj < 0.05, plus absolute log2FC > 0.5 for mRNA). Process NES from preranked fgsea, BH over 1,145 sets.

**(D)** Percent of colocalized genes whose best gene-study pair places at least 50% of its primary-signal posterior in open chromatin of each lineage, in at least one of two liver snATAC sources (circle) or both (triangle). Direct MASLD/liver-fat genes, n = 22 (28 pairs); liver-enzyme genes, n = 424 (739 pairs); ten genes are in both. Intervals, gene-bootstrap 95%. Diamonds, expectation with posterior spread uniformly over tested variants. Overlap does not assign a target gene.

**(E)** 45 of 117 programs; all 117 are in the source table. Upper block, 31 programs with a member gene open in its own lineage through a panel D pair. Lower block, 14 programs shown only for a robust protein association, so selected on the plotted statistic. Bars: dark, open through a direct MASLD/liver-fat pair (RORA); light, through liver-enzyme pairs only. Right, partial Spearman rho (participant-bootstrap 95% CI) between program protein score and histologic burden in 46 MASLD participants. Robust: Freedman-Lane BH q < 0.05 over 95 scorable programs and stable sign across score variants and both batches. Gray outline, CI excludes 0 but not robust.

**(F)** Left, illustrative H&E overlay of the adjusted Stromal ECM score (section JBO019, donor H37), chosen outcome-blind; no pathologist annotation. Right, Moran excess over matched genes after adjusting spot scores for all 16 cell2location abundances, in Guilliams et al. and Vu et al. [ref: Vu et al., liver spatial transcriptomics]. Small points, donors (Guilliams et al., 4) or arrays (Vu et al., 10); band, matched-null 95% range. Labels give BH q from 9,999 matched-gene draws (two programs per source) and donor or array counts. The gene null holds donors and arrays fixed, so it is not a donor-level test; Stromal ECM excess is positive in 3 of 4 Guilliams et al. donors (exact one-sided sign-flip p = 0.19; the smallest possible p with 4 donors is 0.0625). IGFBP7 carries 2.47% of the Stromal ECM weight and BICC1 5.44% of the Ductular injury weight, so the result concerns program scores, not either gene. Vu et al. has no donor key (open point, source-dependent). The Guilliams et al. donors are surgical samples taken for other conditions, not a MASLD cohort; H37 and H38 are also in the single-cell atlas.

## Figure 6

**Figure 6. The MASLD Gene Catalog keeps eQTL-colocalization and bulk RNA calls separate for 28,233 genes.**

**(A)** Catalog structure. Each gene identity carries a separate eQTL-colocalization call and bulk RNA call, each qualified by testability and provenance, plus one molecular object (gene, RNA or DNA element) and one proposed discriminating experiment. Protein, chromatin and spatial context keep their assay or program unit and enter no class or route.

**(B)** Genetic call (rows) by bulk RNA call (columns) for 28,233 gene identities. The genetic call is supported when SuSiE-COLOC PP.H4 > 0.5 in any of 35 Tier-1/2 GWAS (15 direct MASLD, 20 liver enzyme). The RNA call is supported when BH q < 0.05 (over 24,115 genes) and |log2FC| > 0.5 in one limma-voom fit of 774 genetic-identity clusters (138 control, 636 disease), which are not verified donors. Blue cells are tested in both domains; gray cells are untested in at least one. No negative-call rule was defined, so there is no tested-negative row or column, and untested does not mean negative. Genetics untestable means no Tier-1/2 COLOC record (9,312) or SuSiE-untestable in one study (7; finite PP.H4 elsewhere, maximum 0.164). The untestable-by-untestable cell is 0 by construction. "Coding PIP ≥ 0.5" counts 31 genes with a fine-mapped coding variant at SuSiE PIP ≥ 0.5. It does not enter the genetic call, so PNPLA3, TM6SF2 and MTARC1 stay unresolved.

**(C)** Four worked entries. Columns give SuSiE and ABF PP.H4 with the top study, any fine-mapped coding variant, bulk RNA log2FC with 95% CI (log2FC ± 1.96 SE; 774 genetic-identity clusters) with BH q, and the Catalog class. The line under each gene gives its gene-level and molecular-object routes. Dashed lines mark |log2FC| = 0.5. ABF PP.H4 is shown for reference and makes no call. *GLP1R* is untested: RNA not detected and no Tier-1/2 SuSiE result. *MTARC1* has no Tier-1/2 SuSiE result, so its genetic call is indeterminate. *MTARC1* and *THRB* have indeterminate RNA calls because the change is below the 0.5 floor. The *MTARC1* missense variant (PIP 0.76) does not enter its genetic call. "Both supported" (*HKDC1*) does not compare direction. Genetic support does not establish a causal mechanism.

**(D)** Genes per Catalog class, with the gene-level route printed beside each class: untested 13,551, unresolved 13,127, established-state associated 1,129, genetically anchored 400 and supported on both axes 26 (Catalog class "concordant"). The untested bar splits into no RNA call (4,232, not detected) and no genetic call (9,319). Gray is untested and blue tested. Hatching marks the 708 untested genes with one supported call (698 RNA; 10 genetics, narrower than a line). The 20 unresolved genes with a coding variant at PIP ≥ 0.5 are routed by their eQTL-colocalization class. Routes are proposals, and counts are not a priority ranking.

**(E)** Ways to use the Catalog: view the worked entries, inspect each molecular object with its uncertainty, read source references, and download the gene calls with the molecular-object records. All views draw on the same inputs and give no combined score or rank. The Catalog files are available at [accession to be supplied] (see Data availability).

## Figure 7

**Figure 7. Predictive models are scored against task-specific baselines beside a record of their training exposure.**

**(A)** Evaluation framework. Each of 8 registered tasks has a baseline set in advance. A score whose exposure record is unresolved is uninterpretable, not evidence of skill. The bar counts audit entries by checkpoint finding (104 audits of 106 model entries), before eight later audit revisions. These checkpoint findings are separate from the study-specific exposure in (C) and (F).

**(B)** Gain over each task's own baseline, named with its score on each x axis. Reporter assay (GSE281364, HepG2 contexts only): Bonferroni-corrected intervals across the seven rows (paired block bootstrap); sequence-based reporter baselines were not run. eQTL |z| and sign: held-out LD fold, paired 95% LD-block bootstrap. Fibrosis: samples are the unit because the training pool has no donor key; 95% sample bootstrap. H3K27ac: 95% participant bootstrap. ATAC: donor jackknife.

**(C)** Donor-class-balanced macro F1 for broad cell classes on all seven held-out studies (hollow; 102 donors, 50,000 cells) and on each encoder's clean studies (filled). Gray: HVG+PCA with the same head. Intervals are 95% study-stratified donor-cluster bootstraps. Geneformer V2 has no clean study, so its scores are not evidence of skill. Hollow and filled scores use different study sets and are not compared with each other; (F) compares each encoder with HVG+PCA on the same studies.

**(D)** Signed eQTL: best head's sign AUC per |z| stratum (95% LD-block bootstrap, uncorrected across five strata) beside within-block permutation chance and the sequence-attributable bound, which assumes at most one causal variant per gene and LD block. Severity transfer: released instrument, no refit; Spearman with stage (stage-stratified 95% bootstrap) beside each cohort's stage tie ceiling. GSE213621 units are genetic identities, not verified donors. GSE276114 includes participants with HCC, sampled from peritumoral liver without a per-sample label. LD-block join: Pearson r of each row's effect with its gene's mean effect, own LD block left in or held out. Routing: fitted rule and oracle router against the pre-declared floor.

**(E)** Each statistic is oriented so that positive is the prespecified direction and placed on one scale, (estimate − null)/(threshold − null): 0 is the null, 1 the threshold its verdict uses (minimum detectable effect for severity arms and the zone task; h2 z = 7 for GWAS rows; twice the test's resolution for signed eQTL; the pre-declared floor for routing). Severity-arm intervals (95% stratified bootstrap) show precision only; verdicts use the label permutation null, and an underpowered arm is untestable, not negative. Signed eQTL direction, recorded as tested negative, is drawn untestable because its sequence-attributable bound (open diamond) is below the threshold. Treatment response uses the GSE83452 diet arm (6 resolved vs 16 persisted; permutation p 0.23), which its source analysis records as negative and underpowered.

**(F)** Top: exposure of each held-out study (donors in parentheses) for each encoder. Middle: macro F1 gain over HVG+PCA with the same head, on all seven studies (hollow) and clean studies (filled). Bottom: encoder minus the scVI liver encoder trained here, on the encoder's clean studies. Intervals are 95% donor-cluster bootstraps, uncorrected across rows. In the unweighted study mean, a study with 2 donors and 149 cells gets one seventh of the weight; two of twelve model-by-head blocks change sign when weighted by donors, so no sign claim rests on the unweighted mean.

# STAR Methods

## Key resources table

Versions are those recorded in the analysis environments and run logs.

| Reagent or resource | Source | Identifier |
|---|---|---|
| **Deposited data: bulk liver RNA-seq** | | |
| Five-cohort disease and stage analyses | [ref: source study for each cohort]; GSE130970 controls from Hoang et al. (2019) | GEO: GSE126848, GSE130970, GSE135251, GSE162694, GSE213621 |
| Repeat-biopsy pairs (source-labelled; identity not verified) | [ref: GSE193066 source study] | GEO: GSE193066 |
| Program analyses and program-pair discovery | [ref: GSE174478 source study] | GEO: GSE174478 |
| Program-pair test cohorts | [ref: GSE281797 and GSE225740 source studies] | GEO: GSE281797, GSE225740 |
| Retrospective program transfer | [ref: GSE240729 source study] | GEO: GSE240729 |
| Adolescent liver expression | [ref: GSE66676 and GSE185051 source studies] | GEO: GSE66676, GSE185051 |
| Predictive-model evaluation cohorts | [ref: source studies] | GEO: GSE268273, GSE276114, GSE267145 |
| Postmortem liver expression and histology | GTEx Consortium [ref: GTEx v8] | GTEx v8 |
| Provenance record only (library preparation confounded with case-control status) | [ref: PRJNA512027 source study] | SRA: PRJNA512027 |
| **Deposited data: single-cell, chromatin and spatial** | | |
| Single-cell and single-nucleus liver RNA-seq | Gribben et al., 2024 (GSE202379) [ref]; Andrews et al., 2022 (GSE185477) [ref]; Filliol et al. (GSE174748) [ref]; Guilliams et al., 2022 (GSE192742) [ref]; [ref: other source studies] | GEO: GSE202379, GSE244832, GSE174748, GSE185477, GSE189600, GSE192742, GSE136103 |
| Independent single-nucleus RNA-seq | [ref: source studies] | GEO: GSE256398, GSE289173 |
| Single-nucleus ATAC-seq | [ref: source studies] | GEO: GSE244832, GSE281367 |
| Bulk ATAC | [ref: source study] | GEO: GSE296875 |
| Visium and spatial transcriptomics | Guilliams et al. (GSE192741) [ref]; Vu et al. [ref]; HMSMA [ref]; Yakubovsky et al. (2026) [ref] | GEO: GSE192741, GSE287826; GSA-Human: HRA007511 |
| CosMx panel | Boesch et al., 2026 [ref] | GEO: GSE312698 |
| Massively parallel reporter assay | [ref: GSE281364 source study] | GEO: GSE281364 |
| **Deposited data: proteomics** | | |
| Liver and plasma DIA-MS | [ref: PXD051911 source study] | PRIDE: PXD051911 |
| Plasma proteome (excluded; labels not authenticated) | [ref: PXD052937 source study] | PRIDE: PXD052937 |
| **Deposited data: genetics** | | |
| Million Veteran Program Release 4 GWAS | [ref: MVP Release 4] | dbGaP: phs002453.v1.p1 |
| Published liver-trait GWAS | Source studies | PMIDs 31311600, 32298765, 34128465, 34841290, 34957434, 36402844, 36280732 |
| FinnGen GWAS | [ref: FinnGen] | FinnGen R12 |
| UK Biobank liver-enzyme GWAS | Sinnott-Armstrong et al. (PMID 33462484) | GWAS Catalog: GCST90019492, GCST90019497, GCST90019507 |
| Biobank Japan liver-enzyme GWAS | Kanai et al. (PMID 29403010) | NBDC Human Database: hum0197.v3 (BBJ.ALT.v1, BBJ.AST.v1, BBJ.GGT.v1; GRCh37) |
| Pan-UK Biobank round 2 | [ref: Pan-UK Biobank] | Biomarker codes 30620, 30650, 30730 |
| Liver eQTL meta-analysis | Broadaway et al. (2024) [ref] | Liver_eQTL_Meta_Marginal_20240530 marginal summary statistics (hg19) |
| Liver caQTL | Currin et al. (2025) [ref] | doi:10.1101/gr.279741.124 (PMID 40393811); caQTL summary statistics (hg38) |
| European LD reference | PolyFun UK Biobank LD [ref] | UKBB_LD, broad-alkesgroup-ukbb-ld bucket (Amazon S3) |
| Non-European LD panels and ancestry projection | 1000 Genomes [ref] | Phase 3 integrated release v5a (20130502; AFR, AMR, EAS and SAS subsets; GRCh37) for LD; high-coverage Illumina phased panel (founders) for ancestry projection |
| **Annotation resources** | | |
| Gene annotation and genome | GENCODE [ref] | GENCODE v49 on GRCh38.p14; GENCODE v26 (GTEx); GENCODE v46lift37 basic (hg19) |
| Constraint | gnomAD [ref] | v4.1 |
| Gene essentiality | DepMap [ref] | Public 24Q4 |
| Secretome and blood concentration | Human Protein Atlas [ref] | [version to be supplied] |
| Essential and nonessential gene sets | [ref: CEGv2 and NEGv1] | CEGv2, NEGv1 |
| TF regulons | DoRothEA [ref] | Confidence A/B/C |
| Gene sets | MSigDB Hallmark, KEGG, Reactome, GO [ref] | MSigDB v2025.1 (msigdbr 25.1.1): H, C2 CP:KEGG_MEDICUS, C2 CP:REACTOME, C5 GO:BP, GO:MF, GO:CC |
| Liver ABC enhancers | [ref: liver ABC-enhancer source] | AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3; hepatocyte-ENCODE, HepG2-Roadmap and liver-ENCODE biosamples, non-promoter elements (GRCh37) |
| Histology-anchored continuum signature | Kamzolas et al. (2026) [ref] | 145 genes |
| Hepatocyte protein zonation | Weiss et al. (2026), Supplementary Table 3 [ref] | n/a |
| Published MASLD gene lists | Kumak et al. (2026) [ref] | n/a |
| **Software and algorithms** | | |
| STAR | [ref: STAR] | 2.7.11b |
| featureCounts | [ref: featureCounts] | Subread 2.1.1 |
| edgeR | [ref: edgeR] | 4.4.2 |
| limma (voom, voomWithQualityWeights, treat) | [ref: limma]; McCarthy and Smyth (2009) | 3.62.2 |
| ashr | [ref: ashr] | 2.2-63 |
| nlme | [ref: nlme] | 3.1-168 |
| susieR | [ref: SuSiE] | 0.14.2 |
| coloc | [ref: coloc and coloc-SuSiE] | 5.2.3 |
| PolyFun | [ref: PolyFun] | [version to be supplied] |
| OTTERS | [ref: OTTERS] | GitHub daiqile96/OTTERS, commit 26098ad |
| Ensembl VEP | [ref: VEP] | [version to be supplied] |
| liftOver | [ref: liftOver] | rtracklayer 1.66.0 with the UCSC hg19ToHg38 chain |
| HGNChelper | [ref: HGNChelper] | 0.8.15 |
| Hotspot | [ref: Hotspot] | hotspotsc 1.1.3 |
| scVI | [ref: scVI] | scvi-tools 1.4.2 |
| decontX | [ref: decontX] | celda 1.26.0 |
| Scrublet | [ref: Scrublet] | scanpy 1.12 (`sc.pp.scrublet`) |
| LIANA | [ref: LIANA] | 1.7.1 |
| decoupleR | [ref: decoupleR] | 2.12.0 (R) |
| BayesPrism and InstaPrism | [ref: BayesPrism]; [ref: InstaPrism] | InstaPrism 0.1.6; BayesPrism [version to be supplied] |
| cell2location | [ref: cell2location] | 0.1.5 |
| Leiden clustering and UMAP | [ref: Leiden]; [ref: UMAP] | rapids_singlecell 0.14.1 (single-cell atlas) |
| ScaleSC with Harmony | [ref: ScaleSC]; [ref: Harmony] | ScaleSC 0.1.0 |
| CellTypist | [ref: CellTypist] | 1.7.1, Healthy_Human_Liver model |
| PLINK | [ref: PLINK] | 1.9 |
| SnapATAC2 and MACS3 | [ref: SnapATAC2]; [ref: MACS3] | SnapATAC2 2.8.0; MACS3 3.0.4 |
| TranscriptFormer (TF-Sapiens), Geneformer V2, BulkFormer, HyenaDNA, Nucleotide Transformer, UCE, SCimilarity | [ref: each model] | Checkpoints as stated in Method details |
| R and Python | R Project; Python Software Foundation | R 4.4.3 (bulk analyses), R 4.5.3 (decontX); Python 3.10.18 (bulk quantification), Python 3.12.13 (single-cell and spatial analyses) |

## Resource availability

### Lead contact

Requests for further information should be directed to the lead contact, [name and email to be supplied].

### Materials availability

This study did not generate new materials. All analyses used previously generated human data.

### Data availability

This study analyzes publicly available human datasets. Accessions are listed in the key resources table and in the relevant Method details subsections, and GWAS sources are listed in [supplementary table of GWAS sources; number to be assigned]. Single-cell, ATAC and spatial inputs include GSE192742 (Guilliams et al., 2022), GSE244832 (single-nucleus RNA and snATAC), GSE281367 (snATAC), GSE296875, GSE192741, GSE287826 and HRA007511 (HMSMA). The MASLD Gene Catalog gene table (28,233 rows and 126 columns, none of which is a combined score or rank), the 28,233 molecular-object records and the column dictionary are deposited at [accession to be supplied]. The deposit will also include source checksums, access dates, the testing universe of each analysis, a data dictionary and the source-dependency graph used to build the Catalog. No external web portal is published.

### Code availability

Code for every analysis reported in this paper is available at [repository URL to be supplied] and archived at [archive DOI to be supplied]. R analyses used a single recorded conda environment. Python spatial, benchmark and Gene Catalog analyses used their own recorded environments [environment files to be supplied with the archive]. Every randomized analysis records its seed and its complete permutation family. Seeds are stated with each analysis below.

## Experimental model and study participant details

### Human data sources

All analyses reused previously generated human data from public repositories and published summary statistics. No new participants were recruited and no new samples were collected. [Ethics statement for reuse of public human data to be supplied.] Cohort composition, assay and inclusion rules are given with each analysis in Method details.

### Units of inference

Participants or biological donors define replication throughout. Cells, spots, sequencing runs, guides, arrays, sections and repeated measurements never add to the biological sample size. Where a public deposit is keyed on sequencing run, we collapsed runs to donors with the recorded donor-pairing tables before fitting any model. Assay-specific sample totals, including repeated biopsies, were not summed into a total of unique people.

In the bulk RNA-seq cohorts, the unit is a source-qualified genetic-identity cluster with one selected library per cluster. Genotypes called from the RNA-seq reads identified repeated bulk libraries and conflicting sample labels. However, these RNA-based calls do not give complete identity confirmation, and calls in GSE193066 remain uncertain. Genetic-identity clusters are therefore not verified independent donors, and fits on source rosters are conditional reconstructions. A separate identity sensitivity (seed 20260929) excluded ambiguous calls and every group with conflicting labels, and kept the lexicographically first eligible library in each concordant same-person component. It changed no donor roster or gene membership. Inference that depends on these units requires the identity evidence to be resolved, because resampling cannot correct an incorrect biological unit.

Expression-based sex and PCA checks flag discrepancies. They do not establish genetic sample identity or justify relabeling participants. We did not perform a complete cross-study check for overlapping participants.

## Method details

### Bulk RNA-seq cohorts and quantification

The bulk analysis combines five control-bearing cohorts (GSE126848, GSE130970, GSE135251, GSE162694 and GSE213621). It uses 774 source-qualified genetic-identity clusters with one selected library each: 138 controls and 636 disease units. These clusters are not verified independent donors (see Units of inference). Paired-end fragments were aligned with STAR and quantified at gene level with reverse-stranded featureCounts against GENCODE v49 on GRCh38.p14 [ref: STAR] [ref: featureCounts]. Native fragments were recounted against the resolved gene annotation. The 855 annotations that could not be tested were counted but excluded from testing. Count and metadata columns were joined through the selected roster before any model was fitted. Fragment counts, individual read counts and single-cell UMIs are distinct measurement units and were never mixed. PRJNA512027 was kept as a provenance record only, because its library preparation is confounded with case-control status.

PCA and library-size quality control defined technical inclusion. Samples that failed only the expression-based sex check stayed eligible, and every bulk model adjusts for expression-inferred sex rather than the reported sex field. Under this rule, six disagreements between reported and inferred sex remained in the pooled model. [Supplementary denominator table; number to be assigned] reports cohort and sample keys, count-column alignment, overlapping QC flags and the losses at each analysis step. Denominators from earlier versions of this roster are listed in [supplementary table of earlier roster denominators; number to be assigned].

### Disease-versus-control differential expression

We fitted the primary model in the 774 genetic-identity clusters with limma-voom quality weights and the design `~ dataset + inferred_sex + group_binary`, with dataset as a fixed effect [ref: limma] [ref: voomWithQualityWeights]. Genes were filtered with `edgeR::filterByExpr` under `~ 0 + group_binary`, counts were RLE-normalized, and weights were estimated with `limma::voomWithQualityWeights` [ref: edgeR]. The complete tested family contains 24,115 Ensembl rows. A gene was called disease-associated when its Benjamini-Hochberg (BH) adjusted p was below 0.05 and |log2FC| exceeded 0.50. This rule gave 1,854 calls (1,333 higher and 521 lower in disease). As a comparator on the same fit, the effect-size-aware interval-null FDR gate (`limma::treat` at lfc = 0.25, BH FDR < 0.05) gave 2,131 calls. Cells and repeated libraries never enter as independent units.

We report both thresholds because a fold-change filter applied after a point-null test does not in general preserve FDR control within the filtered subset. It also does not establish that the true effect exceeds the chosen floor. The interval-null test builds its effect threshold into the test (McCarthy and Smyth, 2009) [ref: McCarthy and Smyth, 2009]. Its 0.25 threshold is distinct from the 0.50 threshold used for calls and display. Duplicate gene symbols were resolved without reference to outcome. Raw point-null results and ashr shrinkage estimates are sensitivity analyses only [ref: ashr].

GSE130970 healthy controls are the six biopsies reconstructed from Hoang et al. (2019) Supplementary Table 1, not every biopsy at fibrosis stage F0 [ref: Hoang et al., 2019]. The 774-unit model keeps these control labels. Results from fits with earlier control definitions are given in [supplementary table of fits with earlier control definitions; number to be assigned].

### Per-cohort disease models and GSE126848

Per-cohort disease models kept the pooled 24,115-gene universe and normalization, re-estimated voom quality weights within each cohort and used `~ inferred_sex + group_binary` when sex varied. A constant sex term was omitted. Missing sex, or a disease contrast that could not be estimated, was recorded as such. GSE126848 was modeled as `~ group_binary` in 30 all-male genetic clusters (23 control, 7 disease). We reproduced its 6,147 calls independently with a separate group-means basis and checked them against the native selected fragment counts. The identity-site screen removed 24 of 31 disease records, including all ten women in the source, and 2 of 25 controls. This unequal eligibility may limit how well the retained units represent the source cohort, and an estimable design does not resolve that. Retained disease libraries are shallower than controls. Figure 3C shows the unadjusted `~ group_binary` calls. Three sensitivities are reported in the text (Supplementary Note 7) and are not used for the panel: refits of `~ group_binary + log2(library size)` with library size counted before or after RLE normalization; the number of calls on genes with zero reads in all 30 libraries; and a refit after `filterByExpr` within the 30 units with RLE factors recomputed (3,440 calls on 15,592 genes, a different BH family).

### Postmortem liver fat in GTEx

In GTEx v8 liver [ref: GTEx v8], fat was measured on each donor's H&E slide as the clear-vacuole area fraction of the tissue. Thresholds were set before any expression result was read. The gene family, set before the run, was the 1,002 disease-versus-control calls present in GTEx (GENCODE v26) that pass edgeR `filterByExpr` on GTEx liver counts. For each gene, limma `voomWithQualityWeights` (weights from the covariate-only design) and `lmFit` estimated the slope on fat rank in 222 donors with `~ fat_rank + age_group + sex + Hardy scale + ischemic time + RIN + collection center`. The primary statistic was the Spearman correlation across genes between this slope and the bulk disease-versus-control log2 fold change. It was tested two-sided by Freedman-Lane permutation within six strata of center by Hardy group (10,000 draws) [ref: Freedman-Lane permutation]. Two qualifiers formed a Holm family of two: the same correlation in the 138 donors with Hardy scale 1 to 4 (deaths not on a ventilator), and the difference between the fat and ischemic-time correlations in a joint model (2,000 bootstrap resamples). Positive controls were PLIN2 and 197 steatosis-associated genes from the bulk cohorts that are not disease-versus-control calls. Negative controls were slide tissue area and ischemic time. Donors are the unit of replication. The seed was 20260923.

### Adolescent liver cohorts

GSE66676 provides liver expression (Affymetrix HuGene 1.0 ST) from 67 adolescents with severe obesity, sampled at bariatric surgery [ref: GSE66676 source study]. The signature was the 1,854 disease-versus-control calls, of which 881 up and 183 down gene symbols map to array probes. Each probe was z-scored against the No NAFLD group, z was averaged over the probes of each symbol, and the score was the mean z of up symbols minus the mean z of down symbols. We removed two participants whose recorded sex did not match their expression-inferred sex, leaving 32 with NAFLD and 33 without. The single prespecified test (no correction) was the NAFLD coefficient in `score ~ nafld + sex + age + bmi` (residual df 60). We checked it by Freedman-Lane permutation (10,000 draws) and report it as an adjusted standardized difference d with its 95% interval. We computed the minimum detectable d at 80% power for the analyzed n. Random gene sets matched on No NAFLD SD decile served as negative controls. GSE185051 (52 children) was used only descriptively, by correlating the score with NAS [ref: GSE185051 source study]. Participants are the unit of replication. The seed was 20260923.

### Composition sensitivity of disease-associated genes

We refitted the 774-unit disease-versus-control model with estimated cell-type fractions added as log-ratio coordinates. We used four prespecified codings (16 cell types with a floor of 1e-6 or 1e-4; seven lineage groups; and a five-group coding in the three cohorts where it resolves stroma, 537 units) and two fraction inputs (BayesPrism fractions estimated from the fragment counts used here and from an earlier count matrix of the same cohorts), giving eight combinations. Ratios of adjusted to unadjusted log2 fold change were summarized with 1,000 bootstrap draws that resampled units within cohort-by-group strata (seed 20260923). Per-gene labels required four checks to pass. Both positive controls had to pass (12 Stromal ECM calls and 159 stromal-specific calls, each with a median ratio below 0.5 in every combination). The negative control had to pass (52 calls expected to change uniformly, median ratio at least 0.7 in every combination). A within-hepatocyte benchmark and a calibration null (lambda at most 1.2) also had to pass. The stromal-specific positive control did not pass, so no per-gene labels were produced.

### Genetic ancestry sensitivity

Genetic ancestry principal components (PC1 to PC4) were computed from RNA-called genotypes projected onto 1000 Genomes founders, and attached by run identifier to the 774 genetic-identity clusters (138 control, 636 disease) [ref: 1000 Genomes]. We refitted the stored filtered, RLE-normalized DGEList in two ways. T1 used `~ dataset + inferred_sex + PC1 + PC2 + PC3 + PC4 + group_binary`. T2 used the primary design in the 614 clusters called EUR (97 control, 517 disease). limma `voomWithQualityWeights`, `lmFit` and `eBayes` were rerun, and calls used BH over the 24,115-gene universe with adjusted p < 0.05 and |log2FC| > 0.5. The statistic was the number of the 1,854 primary calls that changed status (lost plus gained).

The T1 reference was 200 refits with the four PC columns permuted jointly within each of the 10 cohort-by-group cells, which keeps the PC imbalance between groups. The T2 reference was 100 refits on random subsets drawn within cells to the EUR counts. One-sided empirical p values used (1 + k)/(B + 1), and the family F1 = {T1, T2} was Holm-corrected at family-wise alpha 0.05. Controls were reproduction of the primary fit; the rs738409 G allele frequency by ancestry call (exact CMH test stratified by disease group); a joint F-test of the four PCs in a cohort-by-group cell-means fit, compared between true and permuted PCs; and GSTM1 as a gene expected to differ by ancestry. A gene was labeled ancestry-sensitive, descriptively and without multiplicity correction, if it changed status with the PCs and in at most 9 of 200 permuted fits. We estimated the expected number of false labels by applying the same rule to each permuted fit against the other 199. The seed was 20260923 (L'Ecuyer-CMRG), and all permutations and subsets were drawn before any fit.

### Spike-in calibration of bulk evidence states

We added known disease-versus-control changes to the 774-unit fragment counts (79,095 genes) by binomial thinning (Gerard, 2020) [ref: Gerard, BMC Bioinformatics 2020]. Each count in the group to be lowered was replaced by Binomial(y, 2^-|delta|). In each of 50 replicates, disease labels were first permuted within cohort, keeping each cohort's control and disease counts. We then drew 320 genes per AveExpr decile of the primary fit from the 24,115-gene universe. These genes received uniform changes of 0, 0.5, 0.75, 1.0 and 1.5 log2FC; heterogeneous changes of 1.0 and 0.5 log2FC with between-cohort SD 0.306 (the median DerSimonian-Laird tau among the primary calls); or changes of 1.0 and 2.0 log2FC in one cohort only.

Each replicate was refit with the full primary procedure (`filterByExpr`, RLE, `voomWithQualityWeights`, `~ dataset + inferred_sex + group_binary`, BH over the replicate's universe), together with the five leave-one-cohort-out refits and the five per-cohort fits. A supported call was the primary disease-versus-control call. Source dependence was loss of the call, with the pooled sign, in at least one leave-one-cohort-out refit. We also scored an equivalence rule (DerSimonian-Laird random effects, Hartung-Knapp SE, TOST at margin 0.5 with t on 4 df, BH), but it defines no evidence state. Inside each simulated fit, the unit is the genetic-identity cluster. Spiked genes within replicates are calibration draws, not donors. There is no hypothesis test and no family. Each rate carries a 2,000-draw bootstrap interval that resamples the 50 replicates, reported beside a Clopper-Pearson interval. A threshold certifying a low error rate was met only when the unfavorable bootstrap bound was on the correct side. The pooled, leave-one-cohort-out and per-cohort fits were reproduced before any replicate was run. Stop rules covered the null background (median unspiked supported calls above 105) and thinning size. Seeds were 20260923 + r for replicate r (Mersenne-Twister) and 20261923 for the bootstrap.

### lncRNA annotation and disease association

Versioned Ensembl identifiers from the GENCODE v49 annotation used for counting were the identity key for protein-coding genes and lncRNAs. Symbols were display aliases only. The primary lncRNA analysis excluded small RNAs, circRNAs, pseudogenes and transcript-isoform claims. The five pooled cohorts had to share the same reverse-stranded featureCounts setting. We checked strand orientation from deposited STAR summaries and, for GSE162694, from a deterministic recount of ten BAM files with featureCounts modes 0, 1 and 2. The cohort-level criterion required a median reverse-to-same assignment ratio of at least four, at least 90% of checked libraries at four or above, and every library above one. GSE213621 also required that no clinical stage group have more than 10% of libraries below four, and that the largest between-group difference in failure fraction be at most five percentage points. No sample was excluded on the basis of this check.

lncRNA tests share the 24,115-gene pooled family, including the primary and interval-null results; we did not calculate a separate lncRNA-only FDR. Cohort-specific, leave-one-cohort-out and equal-library-weight fits were prespecified robustness analyses. The named-example and program-removal analyses below were run on an earlier 23,702-gene version of the family and have not been rerun on the current fit. A named lncRNA example also required a canonical chromosome, an unambiguous genomic mapping, concordant direction in at least four cohorts, no materially opposite confidence interval in the fifth cohort, concordant direction in all five leave-one-cohort-out fits, retained interval-null significance in at least four of those fits, and concordant equal-weight results.

The secondary association between biotype and disease-association status used a gene-level logistic model adjusted for average expression, expression variability, log gene length and the number of cohorts in which the gene was detected. Genes are comparison units in this descriptive analysis, not biological replicates. For the direction comparison, each disease-associated lncRNA was matched 1:1 without replacement to a protein-coding gene on these covariates. The caliper was chosen without reference to outcome as the first of 0.20, 0.15, 0.10 and 0.05 SD at which all standardized differences met the 0.10 limit. In the 774-unit analysis, the 0.20-SD match failed this balance rule, so 0.15 SD was selected, giving 520 pairs. Gene-level matching is descriptive.

Each lncRNA was assigned one deterministic genomic-context class from exon overlap, gene-body overlap, strand and opposite-strand TSS proximity to protein-coding genes, applied in that order of priority. Same-strand exonic-overlap and complex mappings stayed in coverage tables but could not become named examples. Genomic context was not interpreted as evidence that the mature RNA molecule is functional.

Prespecified Hotspot programs with at least 5% lncRNA L1 weight were rescored after removal of all uniquely mapped lncRNA members. Module membership of other genes was not changed, and modules were not rediscovered. The exact original cell scores had been saved, but the original neighbor graph had not. We therefore reconstructed the source graph from the original cell-type matrix, source filters, Hotspot model, neighbor count and gene filter. The exact stored score was kept as the baseline, and only the reconstructed gene-removal difference was added to it. The technical reconstruction required a correlation of at least 0.995 with the stored all-gene score. This threshold was set during source checking and kept separate from biological acceptance. Scores were collapsed to biological donor. A sensitivity passed when donor-level Pearson and Spearman correlations were both at least 0.90 and the coefficient direction in the prespecified model `score ~ stage_ordinal + dataset` was preserved.

### Fibrosis-stage and NAS contrasts

We compared adjacent fibrosis stages and adjacent NAS groups cross-sectionally in the five-cohort fragment-count data. Fibrosis contrasts used the three cohorts with an exact F0–F4 fibrosis stage (GSE130970, GSE135251 and GSE162694), with each contrast drawn from 387 staged units. GSE213621 deposits only F0–F1, F2 and F3–F4 bins and was not used in these contrasts. Each analysis kept one library per selected genetic-identity cluster and only the cohorts that carried both compared groups. The design was `~ dataset + inferred_sex + group` when the multi-cohort terms were estimable. Counts were normalized with edgeR and analyzed with `limma::voomWithQualityWeights`. Empirical-Bayes p values were BH-adjusted separately across the complete 24,115-gene family for every contrast. No stage-specific fold-change threshold was applied. Direct F1, F2, F3 and F4 versus F0 models used the same data and covariates. For every model we saved sample count, cohort count, design rank, coefficient name and the complete family, and the five-cohort models passed checks of design rank, family size and FDR recomputation. Arrows in compact panel labels denote ordered cross-sectional contrasts between different participants, who were not sampled repeatedly over time. Raw counts of differentially expressed genes are not treated as comparable effect sizes when sample size or cohort availability differs among contrasts.

The F0 reference combines 14 source controls and 89 disease biopsies (103 clusters). The adjacent F0→F1, F1→F2, F2→F3 and F3→F4 contrasts contain 205, 190, 157 and 94 clusters, and no cluster is on both sides of a contrast. We report this reference with two sensitivities: one removes controls from F0 only, and one adds a control-status term at unchanged sample size. The disease-only-F0 rule does not remove controls at other stages. Contrasts that could not be estimated are recorded as such. Early-stage GNMT and MAT1A claims that do not survive these comparisons are not retained. As a further sensitivity, every fibrosis contrast was refitted with each cohort left out in turn under the same design, and F1 versus F0 was fitted in each cohort alone (`~ inferred_sex + group`). Cohort heterogeneity of F1 and F4 versus F0 was tested with a cohort-by-stage interaction (`~ inferred_sex + dataset + dataset:group`), using a 2-df moderated F test with BH correction over 24,115 genes. Pointwise 95% confidence intervals use the moderated coefficient standard error and gene-specific total degrees of freedom. They are not simultaneous intervals across genes.

### Grouped-stage comparison in GSE213621

GSE213621 records fibrosis only in source bins (F0–F1, F2, F3–F4). To compare it with the exact-stage cohorts, we refitted GSE130970, GSE135251 and GSE162694 on the same bins in disease units only (371 genetic-identity clusters: 65, 197 and 109). This fit used limma `voomWithQualityWeights` on the stored counts and RLE factors with the design `~ dataset + inferred_sex + bin`. GSE213621 was fitted in its 258 disease units (F0–F1 87, F2 94, F3–F4 77) with `~ inferred_sex + bin`, and its 68 controls were excluded. For each contrast (F3–F4 and F2, each against F0–F1), the statistic b was the ordinary least-squares slope of the GSE213621 log2 fold change on the discovery log2 fold change, over discovery BH calls kept in GSE213621. The null permuted bin labels 1,000 times within eight strata of GEO accession quartile by inferred sex, so it keeps submission order and sex structure. P values are one-sided, and the two contrasts form one Holm family (family-wise alpha 0.05). Intervals come from 200 bootstrap draws stratified by bin and sex. The normal interval is primary, and the percentile interval is reported beside it. A contrast was untestable if the minimum detectable slope under the null (2.80 times the null SD) exceeded 0.30, or if fewer than three of the retained positive controls (COL1A1, EPCAM, GNMT, LUM, CYP2C19) had the predicted sign at nominal p < 0.05.

We read effect size from a ratio. Each discovery cohort was held out in turn and the other two were refitted. The mean GSE213621 slope on those refitted effects was then divided by the mean slope of the held-out cohort on the same effects and genes. Similar size required the 95% interval of the ratio to lie within 0.67 to 1.5. The negative control was the slope of the GSE213621 male-versus-female coefficient on the discovery effects over autosomal genes. It fails if its interval excludes 0 and its absolute value exceeds 0.10. Outside the test family, sensitivities added GEO accession quartile, RLE factors recomputed after filtering, a within-cohort top-50-gene share covariate, omission of GSE162694, a women-only fit (161 units), and six isometric log-ratio coordinates of seven lineage groups from the fragment-count deconvolution fractions added to both fits. Per-gene states used BH within the retained genes of each passing contrast. The seed was 20260923. In this analysis, the minimum-detectable-slope rule made F2 untestable, and the negative control failed. The GSE213621 F3–F4 versus F0–F1 contrast had been fitted, and its sign agreement read, before this analysis was specified, so the analysis is not blind.

### Library properties of the stage contrasts

The fibrosis and NAS stage contrasts were drawn from 387 staged fibrosis units and 395 staged NAS units in GSE130970, GSE135251 and GSE162694. For each contrast, we refitted the stage model with three within-cohort standardized covariates: log2 library size, the Gini coefficient of the count distribution and the stored normalization factor. We compared the share of reported BH calls that stayed BH-significant with 200 refits in which the three covariates were permuted together within cohort. F1 versus F0 was a family of one; the other ten contrasts formed one BH family. The seed was 20260923. As a positive control, a planted covariate correlated with stage group (about 0.6) replaced the three covariates for F1 versus F0 in 20 draws, and the retained share had to fall below the null 2.5th percentile in every draw. It did not, so the prespecified stop rule applied and no result from this check is reported.

### Bulk NMF axes

We interpreted bulk NMF programs as continuous axes, not discrete subtypes of participants. Cross-seed matched-factor cosine similarity was used to assess the stability of the factor subspace. Hard sample partitions were not treated as reproducible subtypes, because the sample-level adjusted Rand index was modest and no internal metric uniquely selected six factors.

### Histology-anchored molecular continuum

We adapted the molecular continuum of Kamzolas et al. [ref: Kamzolas et al., 2026; doi:10.1038/s42255-026-01543-7] as an additional evidence layer. It does not replace recorded histology or the disease-versus-control estimand. GSE126848, GSE130970 and GSE135251 overlap the source study's discovery cohorts. GSE162694 and GSE213621 were evaluated separately as cohorts independent of that study. A within-person sensitivity used the 54 of the 58 source-labelled GSE193066 repeat-biopsy pairs in which both libraries passed technical quality control. Their identity has not been verified, so they do not provide primary replication evidence. Continuum and molecular scores were calculated in all five cohorts. Descriptive fixed-window figures show all five cohorts. Stage-stratified figures show GSE130970, GSE135251, GSE162694 and GSE213621, because GSE126848 lacks recorded fibrosis. The source-overlap cohorts were not added to the primary stage- and sex-adjusted inference or meta-analysis. Versioned Ensembl identifiers were reduced to GENCODE v49 base IDs. Of the 145 signature genes, 139 were observed; the six unavailable genes are listed in the mapping table.

The two co-primary scores were a cohort-specific PC1 of the 145 genes and an external projection with prespecified discovery loadings. PC1 signs were oriented to the prespecified discovery direction, not to target histology. A mean percentile rank of the two scores was used for display only. Within GSE162694 and GSE213621, each molecular outcome was modeled separately against each score as `outcome_z ~ continuum_z + factor(fibrosis_stage) + inferred_sex`. Cohort coefficients were combined by fixed-effect inverse-variance meta-analysis. Random-effects estimates and NAS substitution in GSE162694 were sensitivities. Natural splines with three degrees of freedom were compared with the linear model under complete-family correction. Nine prespecified 20-percentile-wide windows with 50% overlap were calculated within each displayed cohort after within-cohort feature standardization. The windows were descriptive and were never used to select features or define significance.

Gene-level models used count-based `limma::voomWithQualityWeights` and the complete 24,115-gene BH family. Of the 145 signature genes, 143 were observed in this 24,115-gene family. These 143 genes were excluded from ordinary self-association inference. They were tested only with target-specific leave-one-gene-out versions of both axes, in a separate 143-gene family. We describe these results as continuum-associated gene expression. All signature genes were removed before pathway, TF-regulon, Hotspot or fixed-basis NMF scoring. A pathway or program was testable only if at least ten genes and at least 80% of its original membership or signed L1 weight remained. Multiple-testing correction used the full prespecified family for each declared score and model: 50 Hallmark, 658 KEGG, 1,787 Reactome, 7,583 GO biological-process, 1,855 GO molecular-function, 1,042 GO cellular-component, 268 TF, 117 Hotspot and ten NMF tests [ref: MSigDB and gene-set collections].

The ten NMF loadings kept their original k4/k6 definitions and continuous interpretation. Inferential NMF results used signature-excluded fixed-W scores only when three conditions held: full-basis reconstruction correlated at least 0.99 with the released H, at least 80% of factor weight remained, and excluded and native scores correlated at least 0.90 in both evaluation cohorts. The 117 Hotspot programs kept their fixed memberships and weights. Paired GSE193066 models used `delta_outcome_z ~ delta_continuum_z + delta_fibrosis`. Because the pair identities are not verified, these models add no evidence. The bulk continuum analyses use selected genetic-identity clusters, one library per cluster, not verified independent participants. We describe results as cross-sectional stage-associated remodeling, not as ordering in time within a person.

### Continuum score sensitivities

We assessed both continuum scores with categorical histology and inferred sex, with and without 15 isometric log-ratio composition coordinates, in three rosters: the original roster, a disease-only roster and the exact-stage GSE162694 roster. The exact-stage roster was also deliberately coarsened, with 2,000 paired within-cohort resamples to estimate coefficient differences. Metadata kept source control status, original histology category, documented exact stage and stage resolution. GSE213621 categories are Control, F0–F1, F2 and F3–F4, and a combined category never supplies an exact stage. We checked categorical relabeling and reproduction of stored coefficients directly. The 117 programs and the 43 molecular systems remained separate BH families, with ordinary and HC3 results kept separate. Age was a complete-case sensitivity, using source-documented age, against a matched baseline. Unavailable BMI and other clinical covariates were not imputed. The fixed 117-program weights and the complete testing families were unchanged.

### Molecular-systems map

For the molecular-systems display, we first kept a complete 13,092-row testability table spanning Hallmark, KEGG, Reactome, three GO collections and the 117 fixed Hotspot programs. A feature was eligible for the map if it passed its prespecified gene-count and retained-membership criteria in both evaluation cohorts. The 139 observed genes from the 145-gene continuum signature were removed before member genes were mapped to GENCODE v49 identifiers. Weighted feature-by-gene membership was inverse-document-frequency transformed, reduced to 100 singular-vector coordinates and used to build a mutual 30-nearest-neighbor cosine graph. Leiden clustering was restricted to the largest graph component [ref: Leiden]. Three disconnected two-node components were recorded as microcomponents and excluded from system summaries. Systems received size-ranked identifiers, not biological names. A seeded UMAP of the membership coordinates was used only for display [ref: UMAP]. Coordinates, component assignments and system memberships were checksummed before any stage, continuum or donor-score outcome was read.

The molecular-systems analysis, including the cluster-level 43-system models behind Figures 4G and 4H, uses the 774-unit, 24,115-gene expression data. In the five-cohort sample table, the GSE213621 stage code places its 68 controls, which have no recorded histology, at stage 0. Figure 4G and 4H stage and continuum coefficients are therefore estimated in disease units only (109 GSE162694 and 258 GSE213621). The system maxT and competitive-null tests and the supplementary comparisons of stage and continuum layers were run on 435 units. The per-layer effect correlations were also refitted in the 367 disease units, without recomputing the layer counts. In the stage factor of the continuum model and in the maxT strata, those controls form their own stage category. In the linear stage-trend axis, they enter as stage code 0. The system redundancy summary used a 466-unit descriptive roster. Fixed system memberships and outcome-blind coordinates were unchanged.

Fibrosis-stage and stage-adjusted fixed-projection coefficients were overlaid on the same coordinate table. System summaries were calculated in two ways: as the equal mean of within-collection medians, and as the median across all member features. Fixed-window summaries used the nine prespecified windows in all five cohorts under both aggregation rules. Cross-system correlations, direction counts and window-direction concordance were treated as descriptive, because gene sets are nested, participants recur across molecular features and adjacent windows overlap. We used no independent-sign null, system-level p value, outcome-adaptive community or single-medoid biological system name.

For the pooled molecular-systems window display, we first calculated the collection-balanced system score separately in every cohort and window: the median feature score within each represented collection, followed by an equal-weight mean across collections. We then averaged the five cohort-specific system scores with equal cohort weight, so GSE213621 did not receive greater weight for containing more participants. Continuum percentiles and feature standardization remained cohort-specific, so this display is an aligned descriptive summary, not a jointly normalized participant pool. For display in Figures S4T and S4U, systems were sorted by the descriptive equal-cohort pooled late-minus-early contrast: the mean score in windows centered at percentiles 80 and 90 minus the mean score in windows centered at 10 and 20. The five cohort-specific heatmaps use the same row order. This result-dependent ordering was not used for feature selection, support calls, nonlinear timing or any inferential test.

### Continuum-linked direction, programs and readouts

This extension asks two questions. Do inherited expression associations and tissue remodeling point to different experiments? Which measurements could describe molecular variation within an existing histologic eligibility group? The analysis is exploratory. It does not change genetic memberships or Gene Catalog evidence classes. It makes no recommendation for clinical practice, does not estimate treatment benefit and does not set a clinical cutoff. The continuum and its sliding-window adaptation are those of Kamzolas et al. [ref: Kamzolas et al., 2026].

*Signed inherited direction.* The signal-level direction analysis used the same 50-study colocalization results as the gene-level analysis. We kept all 5,362 exported SuSiE signal pairs and all 22,356,697 variant posterior records from 816 main-trait gene–study pairs. Signal support used the same strict PP.H4 > 0.5 threshold as the colocalization analysis. Direct MASLD, liver fat, liver enzymes and other traits were kept as separate fields, as were ancestry and the GWAS and eQTL signal identifiers. Supported gene–trait records elsewhere in the portfolio, including ABF-only sensitivity support, stayed visible when signal posteriors were unavailable. They did not take a direction from an unrelated lead variant. Rows without a gene identifier were exported separately and could not be joined to a named Gene Catalog record. Main-trait and full-portfolio testing denominators were recorded separately, distinguishing missing SuSiE tests from tests without support. These posterior support thresholds are not frequentist p values.

At matched GRCh37 chromosome and position, we matched the eQTL effect and other allele pair to the recorded GWAS effect and other allele pair. An exact swap negated the eQTL beta. Unmatched alleles, non-SNV alleles, ambiguous duplicate effects, missing effects and palindromic variants were left direction-unresolved. We did not infer strand from an effect sign. This palindrome rule affects direction only and does not remove variants from the colocalization analysis. For each shared-variant posterior vector, the sign of GWAS beta × aligned eQTL beta describes expression associated with a higher value of the named trait. Only for direct disease phenotypes do we describe this as expression associated with disease risk. We summed the original shared posterior mass carrying each sign and required at least 0.95 of the original total in one resolved direction. Missing or ambiguous mass was not removed by renormalization. This is a declared descriptive direction criterion, not a calibrated posterior probability that expression mediates the trait.

The available signed effects are marginal summary-statistic effects and do not identify the conditional effects of individual SuSiE signals. We therefore resolved directional perturbation records only for loci with one exported GWAS signal and one exported eQTL signal. Records with multiple exported signals keep marginal-direction evidence but remain signal-direction unresolved. Even a single exported pair does not exclude additional undetected signals. No instrumental-variable effect or drug effect was estimated. Per-study GWAS sample sizes and per-variant eQTL sample sizes accompany the signed evidence. We also tabulated consistency among resolved records for the same gene, named trait and ancestry across studies. Conflicting signs remain in their individual records and suspend a common directional perturbation hypothesis. Consistent signs across studies are not called independent replication, because participant overlap between studies is unresolved. Different signal locations within a gene remain distinct possible mechanisms, not a pooled estimate.

*Expression relationship.* This extension was run on an earlier 23,370-gene version of the bulk testing family, not the 24,115-gene family of the primary fit. We joined inherited evidence to that complete family and kept genetics-only genes outside that measured universe. A matching or opposing direction required the two-axis, two-cohort continuum-support rule. Genes with inherited support but no supported remodeling were kept, and genes not observable in bulk were distinguished from tested genes. The 139 observed axis genes keep their separate leave-one-out results and 139-gene BH family. They do not acquire ordinary continuum support from the axis they define. Opposing directions motivate experiments that test both directions and alternative explanations; they do not show compensation. Lifelong genetic associations do not establish the sign or timing of a pharmacological effect. Gene Catalog annotation tables keep separate fields for the proposed transcript, the perturbation hypothesis, cell-context links, the disease-associated profile, the proposed readout and the unresolved experiment. Shared program membership alone does not establish regulation between a target and a different proposed marker. No composite evidence score or target ranking was introduced.

*Within-stage programs and composition.* We kept all 117 fixed program memberships and weights and reused the signature-excluded scores from the selected genetic-identity clusters. The primary model was `program_score_z ~ continuum_z + factor(fibrosis_stage) + factor(inferred_sex)`, fitted separately in GSE162694 and GSE213621 for each of the two co-primary axes. These analyses used 109 and 326 complete genetic-identity clusters, respectively (435 total), and every program has an explicit denominator and failure state. These cohorts had been examined before, so they provide replication evidence but not untouched prospective validation. The three source-overlap cohorts contribute descriptive profiles only.

The composition sensitivity added 15 orthonormal log-ratio coordinates from the 16 estimated BayesPrism lineage fractions [ref: BayesPrism]. Fractions below 1e-6 were floored before taking logarithms, and a normalized Helmert basis supplied the coordinates. We report the frequency of floored fractions by lineage and cohort. Both models used the same complete genetic-identity clusters and outcome units. Axis scores were standardized within the complete cohort, and singular designs were marked untestable. We kept ordinary least-squares coefficients, standard errors and two-sided t tests, and exported HC3 standard errors and t-reference sensitivity tests separately [ref: HC3 estimator]. Fixed-effect inverse-variance summaries required both evaluation cohorts. BH correction covered all 117 programs separately in each cohort, axis and model family, and in each meta-analysis axis and model family, with untestable members counted in the family size. Support required both co-primary meta-analysis BH q values below 0.05 and consistent signs in all four cohort-by-axis fits. A separate, stricter field required BH q below 0.05 in each cohort and each axis. Neither criterion implies 117 independent biological processes. Composition-adjusted associations can still reflect imperfect deconvolution or shared expression measurements, and attenuation does not estimate how much cell composition contributes. We added no nonlinear test; the three-degree-of-freedom spline tests and their complete-family corrections are unchanged. We make no claim about program coherence from this analysis, because the reconstruction and calibration checks that such a claim requires were unresolved.

*Cell and tissue context.* Single-cell source-lineage labels are kept as discovery context. They do not show the cellular origin of a transcript after ambient RNA correction, and they do not independently confirm a bulk RNA result. Signal-level ATAC overlap and program-level ATAC association keep their own lineage, coverage, donor denominators and evidence states. Spatial assay coverage is reported for the complete program family. Formal matched-background tests are confined to the programs originally tested and their stated multiplicity family. Observable programs that were not tested do not receive a negative spatial result, and a program-level spatial result is never transferred to each member gene. Paired-biopsy program and gene analyses keep their within-person unit, complete testing families and measured changes. Repeat biopsies are not treated as technical replicates. This analysis does not establish within-stage, cell-specific genetic regulation.

*Protein readouts.* Adjusted liver DIA-MS disease-association results keep their log2 abundance contrast and full-proteome BH correction. They come from 58 source samples (12 controls and 46 with MASLD), adjusted for acquisition batch, age, BMI and sex. Earlier program-protein tests keep their original subset and identifiers and are reported separately; they are not mapped to all 117 program-level tests. We export complete member-protein observability without treating individual member-protein associations as a test of an entire RNA program.

*Display and reproducibility.* The integrated figure uses the nine overlapping descriptive windows. Its complete program map averages the five cohort window means equally; per-cohort means and participant counts are available. Window membership did not select genes, transitions, thresholds or significance. Native assay effects and their distinct error models remain separate, and repeated RNA-derived features do not count as independent modalities. All scripts set seed 20260906 and record R session information. Checks of allele swaps, ambiguity, posterior mass, complete families, participant denominators and reproduction of source effects accompany the output.

### GWAS portfolio and phenotype tiers

The primary genetic portfolio contained 35 GWAS strata ([supplementary table of GWAS sources; number to be assigned]), all placed in the main portfolio at tier 1 or 2:

- 12 from the Million Veteran Program Release 4 summary statistics by genetically inferred ancestry (open-access dbGaP phs002453.v1.p1; NAFLD, ALT and AST in AFR, AMR, EAS and EUR) [ref: MVP Release 4];
- 9 from published GWAS (PMIDs 31311600, 32298765, 34128465, 34841290, 34957434, 36402844 and 36280732, the last contributing deCODE, Intermountain and UK Biobank strata);
- 2 from FinnGen R12 (NAFLD, NASH) [ref: FinnGen];
- 3 UK Biobank European liver-enzyme strata (Sinnott-Armstrong et al., PMID 33462484; GWAS Catalog GCST90019492, GCST90019497 and GCST90019507);
- 3 Biobank Japan liver-enzyme strata (Kanai et al., PMID 29403010); and
- 6 Pan-UK Biobank round 2 AFR and CSA liver-enzyme strata (biomarker codes 30620, 30650 and 30730) [ref: Pan-UK Biobank].

Tier 1 contained direct NAFLD, NASH and PDFF traits, and tier 2 contained ALT, AST and GGT traits. Tier-3 mixed-etiology sequelae and tier-4 distal nonspecific traits were supplementary. Direct-disease and enzyme results were tabulated separately throughout. The 35 primary strata covered European, African, East Asian, Admixed American and South Asian GWAS populations; the ancestry label refers to the GWAS stratum. The Broadaway liver-eQTL panel is predominantly European [ref: Broadaway et al., 2024], so we describe analyses that pair non-European GWAS with this panel as cross-population rather than fully ancestry-matched colocalization. An expanded 50-GWAS portfolio was analyzed as a sensitivity and is not substituted into main claims.

### GWAS fine-mapping

GWAS SuSiE fits used ancestry-specific LD [ref: SuSiE]. European analyses used the full PolyFun UK Biobank European LD reference with L = 10 [ref: PolyFun]. Non-European analyses used the ancestry-matched GWAS LD panel built from 1000 Genomes where one was available, while the eQTL side remained European. The fine-mapping unit was one converged SuSiE credible set, identified by study, trait, locus, ancestry and credible-set ID. Sets with raw PIP sums outside 0.90 to 1.10 failed the primary reliability criterion. Lead variants had p < 5 × 10⁻⁸ and were clumped with PLINK 1.9 against 1000 Genomes European genotypes (r² 0.1 within 500 kb) [ref: PLINK]. The Biobank Japan strata used distance pruning at 500 kb instead. [Lead-variant procedure for the non-European MVP strata to be supplied.] Each locus spanned the lead position ± 0.5 Mb in all 50 studies.

Fine-mapping used `susieR::susie_rss` on the LD matrix with 1 × 10⁻⁶ added to the diagonal. Quantitative traits used effect sizes and standard errors, with n equal to the total sample size. Binary traits used z scores, with n equal to the effective sample size 4/(1/N_cases + 1/N_controls). The number of effects, iteration limit and purity threshold were the susieR 0.14.2 defaults (L = 10, at most 100 iterations, minimum absolute within-set correlation 0.5). Credible-set coverage was 0.95, and the seed was set to 42 before each fit. Where SuSiE did not converge, PIPs and credible sets were set to missing and the locus contributes no credible set. Credible sets were counted only for main-portfolio tier-1 and tier-2 studies. The GWAS-side SuSiE fit inside the colocalization analysis used L = 10, a minimum absolute within-set correlation of 0.5, at most 500 iterations and residual variance fixed at 1.

### Liver-eQTL colocalization

GWAS strata were integrated with the Broadaway liver-eQTL meta-analysis of 1,183 donors. Single-variant approximate Bayes factor (ABF) colocalization, which uses no LD, was evaluated for every testable gene–GWAS pair [ref: coloc]. Multi-signal SuSiE colocalization was attempted only when four prerequisites were all available: a converged precomputed eQTL SuSiE fit, sufficient overlapping variants, a usable ancestry-specific GWAS LD matrix and a converged GWAS SuSiE fit [ref: coloc-SuSiE]. If any SuSiE prerequisite failed or was structurally absent, the pair received an explicit ABF-fallback or ABF-only state, which we did not interpret as a negative SuSiE result. The prespecified high-resolution support threshold was SuSiE PP.H4 > 0.5. The broader definition was SuSiE or ABF PP.H4 > 0.5. The evidence-class interface uses the SuSiE call. The complete ABF map and the SuSiE-or-ABF union are reported beside it with their own tested denominators. We describe colocalized genes as genetically anchored candidates; colocalization alone does not establish that a gene mediates the trait association.

For every portfolio and trait tier, we counted method availability from finite numeric posterior values, not from method labels. ABF was finite for every row of the colocalization table and spanned 18,976 stable Ensembl genes in both the full portfolio and the primary subset. Availability was as follows.

| Portfolio | Studies | Finite ABF rows | Finite SuSiE rows | Fraction of table rows | Genes with finite SuSiE |
|---|---:|---:|---:|---:|---:|
| Full portfolio | 50 | 948,503 | 21,799 | 2.298% | 4,869 |
| Primary tier 1/2 | 35 | 664,022 | 12,892 | 1.942% | 3,918 |
| Tier 1 direct disease and liver fat | 15 | 284,519 | 771 | 0.271% | 392 |
| Tier 2 liver enzymes | 20 | 379,503 | 12,121 | 3.194% | 3,888 |

Records are unique under the study–chromosome–gene-ID key. There are 600 anonymous records (420 primary), each with an empty gene ID and symbol and an ABF-only state. These account for all 550 additional rows under a collapsed study–gene key. Identified study–gene records are unique. Anonymous records stay in the complete row denominator but are excluded from all gene counts and named-gene counts. Restricting to identified records gives 947,903 finite ABF records in the full portfolio and 663,602 in the primary subset. No anonymous record has a finite SuSiE result.

Finite results include 20,517 negative SuSiE rows carried forward from the earlier run, all checked numerically, and 1,282 rows that passed the shared-posterior qualification of the rerun. The carried negative rows were not rerun, so this denominator does not assert that every negative passed the new shared-posterior check. Among the 1,307 rerun rows, 1,233 remained positive, 49 were eligible at PP.H4 ≤ 0.5 and 25 were untestable. The corresponding primary counts were 767, 38 and 11 among 816 rerun rows. The 25 failed pairs concern 21 genes (the 11 primary pairs concern nine), and all have another finite SuSiE record. After genes that keep a positive call are counted as positive, 17 full-portfolio and seven primary genes keep a priority untestable state. This partial untestability is distinct from a complete absence of finite results. Full-portfolio nonfinite states comprise 921,663 ABF-fallback rows, 5,016 ABF-only rows and the 25 shared-posterior failures; the primary counts are 647,615, 3,504 and 11. Missing or untestable SuSiE results are not negative tests.

At PP.H4 > 0.5, using stable Ensembl identities, the primary portfolio yields 436 SuSiE genes, 822 ABF genes and a 996-gene union. The full 50-study portfolio yields 708 SuSiE genes, 1,181 ABF genes and a 1,494-gene union. Counts from earlier versions of the colocalization tables are reported in [supplementary table of earlier colocalization counts; number to be assigned].

### Accessible-chromatin QTL colocalization

We tested each European tier-1/2 GWAS signal that colocalized with a liver eQTL (SuSiE PP.H4 > 0.5; 581 study-lead signals) against bulk-liver chromatin accessibility QTL (caQTL) summary statistics from 138 donors (118 of European, 17 of African and 3 of admixed genetic similarity) [ref: Currin et al., 2025]. On the GWAS side, we refit the GWAS SuSiE model of the eQTL colocalization (PolyFun UK Biobank European LD, L = 10), reproduced the stored fits and took the log Bayes factor row of the colocalizing signal. For each accessible peak whose center lay within 500 kb of the GWAS lead, per-variant log Bayes factors were computed from the caQTL effect and its variance under a single-effect assumption (prior SD 0.15). GWAS positions were lifted to GRCh38 and matched on position and unordered alleles, and only SNVs with caQTL MAF ≥ 0.01 could match. Variants that could not be matched kept their posterior mass as unshared. The two full rows were compared with `coloc.bf_bf` (coloc 5.2.3; p1 = p2 = 1×10⁻⁴, p12 = 5×10⁻⁶, the priors of the eQTL colocalization). A signal entered when PP.H4 > 0.8 and the GWAS and caQTL lead variants had LD r² ≥ 0.5, the entry rule of Currin et al. A region was linked to a Gene Catalog gene only when, at the same peak, the caQTL row also colocalized with that gene's eQTL SuSiE effect (PP.H4 > 0.8, same priors). For this caQTL-eQTL comparison, variants were matched on position only. eQTL and caQTL donors may overlap, so this three-way support is not independent.

The unit of description is the GWAS signal, and the unit of inference for the comparison is the independent locus (leads merged within 500 kb). One prespecified test (family of one, no correction) compared liver-only loci with platelet-count-only loci by logistic regression. Platelet-count loci came from the MVP European platelet GWAS. Each arm excluded loci with a lead from the other arm within 500 kb. Chromosomes 1 and 19 were excluded because a pilot had examined them. The model was adjusted for −log10 lead p, the number of genome-wide significant variants within 500 kb, the number of testable peaks and −log10 minimum caQTL p. Stop rules required reproduction of the GWAS fits, positive controls at known loci, at least 120 liver-only and 60 platelet-only testable loci, and a rate of 0.10 or below in a negative control that circularly shifts the caQTL values among shared variants. All passed. The negative-control rates were 0.0407 and 0.0513 without the LD entry rule and 0.0233 and 0.0256 with it. Without the LD rule, the prediction for that negative control, 0.03 or below per arm, was missed. The positive controls on chromosomes 1 and 19 had been examined in the pilot and were not blind; the chromosome 5 controls were blind. Locus proportions carry Wilson intervals, and the descriptive lead-level difference used 1,000 locus-cluster bootstrap resamples (seed 20260923). Fits used susieR 0.14.2.

### Joint testing universe and evidence classes

The joint bulk and colocalization universe contains 14,704 stable Ensembl identities with a disease-versus-control bulk test and a tier-1/2 liver-eQTL record. This universe does not require a finite SuSiE result for every gene. Evidence classes were defined only from the disease-versus-control bulk analysis and the primary SuSiE results; proteomic, spatial and single-cell endpoints did not contribute to class assignment. The Figure 1B overlap is computed on the 3,769 genes whose multi-signal colocalization state is positive (SuSiE PP.H4 > 0.5) or tested at PP.H4 ≤ 0.5. The 10,928 genes with only an ABF result and the seven genes with a SuSiE-untestable pair are excluded rather than counted as negatives. We stratified the expected overlap by 10 AveExpr quantile bins (ties broken by Ensembl ID) and tested it with the exact convolution of per-bin hypergeometric distributions. The two-sided p sums the null probability of every total x with |x − E| ≥ |O − E|. The effect measure was the Mantel-Haenszel OR with its Robins-Breslow-Greenland 95% CI [ref: Mantel-Haenszel; Robins-Breslow-Greenland]. The bins, test and effect measure were specified before the outcome was computed. We checked expression control with a prespecified list of 426 non-colocalized genes, matched one-to-one to the colocalized genes on AveExpr by greedy nearest-neighbor matching. All 1,000 prespecified visiting orders returned 13 disease-associated genes against 21.3 expected within the same bins, so the result is labeled as not fully controlled for expression. Power was computed before the outcome from set sizes for an unstratified exact equal-tailed test: ORs between about 0.46 and 1.77 had less than 80% power to be detected (Supplementary Note 3).

The symbol-grouped comparison before identity adjudication contained 14,720 rows. It is a distinct universe, not the denominator of the joint interface. The Gene Catalog reconciles these 14,720 rows to 14,682 jointly eligible identities, with 31 jointly eligible annotation quarantines and seven priority genetically untestable records outside its joint class universe [exact reconciliation in supplementary table; number to be assigned]. Genes absent from an assay universe are labeled untested rather than unresolved. The Catalog keeps one record per adjudicated identity in the union of the bulk fit and the tier-1/2 colocalization tables. Three groups have no record and are not counted in any class: identities quarantined during adjudication (one symbol mapping to more than one Ensembl ID, or an ID absent from GENCODE v49), second Ensembl IDs of a symbol already kept under another ID, and program-member genes outside both assay universes. Six further colocalization IDs absent from GENCODE v49 are not quarantined. Their genetics rows are carried on the record keyed on the same gene's bulk Ensembl ID [listed in supplementary table]. The Catalog therefore holds 1,853 of the 1,854 disease-versus-control RNA calls.

Evidence classes summarize support categories; they are not powered tests of genetic absence. A missing SuSiE result or a partially untestable gene is not evidence of genetic absence, and assay states are kept separately from the support labels. Classes were defined as follows.

- Genetically anchored (labeled Genetics only in Figure 3G): SuSiE PP.H4 > 0.5 and not a disease-versus-control call.
- Established-state associated (labeled Transcriptomics only in Figure 3G): a disease-versus-control call and no SuSiE PP.H4 > 0.5.
- Concordant, supported on both axes (labeled Both in Figure 3G): both SuSiE PP.H4 > 0.5 and a disease-versus-control call. The direction of the genetic effect is not compared with the RNA direction, because the Catalog genetics rows carry no direction.
- Unresolved (labeled Neither in Figure 3G): jointly testable and meeting neither threshold.
- Untested (not shown as a class in Figure 3G): absent from an assay universe, or untestable in one domain. Quarantined identities have no Catalog record and are not counted in any class.

Sensitivity classes substituted the SuSiE-or-ABF union for SuSiE. Genetic trait scope was assigned from the primary SuSiE result as direct disease, enzyme, both or none, and genetic confidence was assigned as SuSiE, ABF-only or none. A broader prioritized gene universe and alternative transcriptomic contrasts were analyzed as sensitivities and are not substituted into main claims.

### Expression-matched comparison of genetic and disease-state classes

This comparison uses the 774-unit bulk RNA fit and the colocalization tables described above. Genetic positivity requires primary SuSiE PP.H4 > 0.5. Disease-state positivity requires BH adjusted p < 0.05 over all 24,115 tested genes and absolute log2FC > 0.5. Exact stable-ID eligibility defines the 14,704-gene joint universe, and the symbol-grouped and annotation-adjudicated Catalog universes are recorded separately. Genetic-only and state-only genes were matched without replacement on bulk average expression with the hardest-first nearest-neighbor procedure and a 0.05-SD caliper. The same complete-case rules yield 281 hepatocyte-specificity pairs and 73 protein-response pairs. The secondary protein outcome is residualized on absolute bulk RNA effect, so that the outcome is not defined through the RNA-based class selection. Genes and matched pairs are comparison units, not independent donors or independent loci.

The original relabel-and-rematch p values and BH q values for the two outcomes are kept only for reproduction. That null removes the observed expression imbalance and is not a calibrated significance test. The paired sign-flip diagnostic preserves the realized matching (100,000 draws, seed 20260907). It is post hoc, reports nominal p values and does not remove dependence among genes or shared genetic loci. The prespecified directional prediction is reported separately from the observed difference, which was in the opposite direction (Supplementary Note 3).

### Signed comparison of predicted expression with the disease contrast

ALT TWAS z per gene came from OTTERS (P0.05 model; Broadaway liver-eQTL weights; UK Biobank ALT GWAS, GWAS Catalog GCST90019492) [ref: OTTERS]. This OTTERS run predates the current colocalization tables and was not regenerated with them. The disease t per gene came from the disease-versus-control limma-voom fit on 774 genetic-identity clusters, which are not verified donors (control/disease: GSE126848 23/7, GSE130970 6/65, GSE135251 10/197, GSE162694 31/109, GSE213621 68/258). The statistic was the gene-weighted mean, over ten AveExpr deciles, of the within-decile Spearman correlation between z and t across 14,703 genes. The null kept genomic order: 10,000 within-chromosome circular shifts of z, with the stratified statistic recomputed on each draw. The 95% CI came from a delete-one-block jackknife over 200 contiguous blocks of genome-ordered genes. The primary test was a family of one at two-sided alpha 0.05 (seed 20260923). A bound was read when the CI lay inside ±0.03 and the shift p was 0.05 or more. A composition qualifier refit the disease contrast with hepatocyte fraction in the design. Five positive and four negative controls, with results predicted before running, passed; none tests the real TWAS-to-disease pairing. The unit of the statistic is the gene, and the shift null and block jackknife carry the dependence among nearby genes. Results are in Supplementary Note 5.

### Comparison with published MASLD gene lists

We matched genes on the NAS- and fibrosis-associated lists of Kumak et al. (2026) to the Gene Catalog [ref: Kumak et al., 2026]. Among the 3,769 genes evaluated by multi-signal colocalization (SuSiE positive or tested at PP.H4 ≤ 0.5), we compared the number of colocalized list genes with its expectation within ten equal-count AveExpr bins from the 774-unit refit. The test used the exact stratified hypergeometric distribution and the Mantel-Haenszel odds ratio with a Robins-Breslow-Greenland interval. The unit is the gene; genes are not biological replicates. The two lists formed one family with Holm adjustment (family-wise alpha 0.05). We ran a 10,000-draw permutation check of the exact p, a reproduction control, a planted-shift detection control, a positive control comparing the lists with disease state, and two negative controls (1,000 expression-matched random lists per Kumak list, and 1,000 within-list permutations of direction labels) (seed 20260923). The positive controls passed, the permutation check agreed with the exact p, and the direction-label control met its prediction. The expression-matched control missed its prediction but stayed inside the stop limits (mean stratified observed/expected 0.931 for NAS against 0.95 to 1.05; fraction with p < 0.05 0.000 for both lists against 0.02 to 0.08). Detection power is given in Supplementary Note 15.

### Fine-mapped DNA architecture and target-transcript biotype

We analyzed fine-mapped DNA and colocalized target transcripts as separate objects. Within credible sets that met the reliability criterion (raw PIP sum 0.90 to 1.10), SuSiE PIP was normalized and divided into mutually exclusive protein-altering, canonical-splice, synonymous-or-UTR, other-noncoding and unresolved masses. Promoter proximity, liver ABC-enhancer overlap [ref: liver ABC-enhancer source] and human liver lineage-accessible chromatin were nonexclusive context annotations and were not summed as probabilities. Promoter and accessible-chromatin context were assigned by GRCh38 position after 1:1 liftOver of the 16,324 sites to be annotated. We did this because lookup by variant key had failed for sites absent from the annotation table, almost all in non-European studies (Supplementary Note 6). ABC-enhancer overlap was assigned by GRCh37 position [ref: liftOver].

For Figure 2C, each COLOC gene contributed one lead variant, the variant with the highest PP.H4 from SuSiE or ABF. Its class came from direct overlap with hg19 UCSC knownGene transcript features (GenomicFeatures), taking the first that applied: coding sequence or splice site (first or last 2 bp of an intron), UTR, promoter (2 kb upstream to 200 bp downstream of a transcript start), intron, then intergenic.

We estimated direct-minus-enzyme differences in coding plus canonical-splice mass within EUR, within AMR and as a fixed-weight EUR+AMR combination (weights 16/26 and 10/26). We used 10,000 shared physical-cluster bootstrap draws over the 81 pooled clusters holding EUR or AMR credible sets (seed 20260923) and 98.333% percentile intervals (Bonferroni over these three members). A member could carry a trait-class call only if it passed three checks. Coverage: in each contributing ancestry, the share of GWAS lead windows with a credible set differed by at most twofold between classes. Power: the minimum detectable difference, 3.236 times the bootstrap standard deviation, was at most the pooled difference (0.195). The 95% interval for the pooled difference drawn in Figure 2D used 10,000 physical-cluster draws (seed 20260817). Calibration: in 200 random splits of East Asian enzyme clusters (East Asian and African for the combined member) into pseudo-classes with the same cluster counts, the 95% interval covered zero in at least 180 splits and the 98.333% interval in at least 190. A difference was called only if, in addition, its 98.333% lower bound was above zero. A member that passed all checks with an interval containing zero was reported as no difference detected. Noncoding-context differences within AMR used the same calibration rule with 99% intervals (Bonferroni over five contexts). Results are in Supplementary Note 4.

The target-transcript analysis used gene-by-GWAS colocalization rows and classified positive targets by versioned GENCODE biotype. ABF and SuSiE denominators were kept separate, as were direct MASLD/PDFF and enzyme-trait results. A credible set was not assigned to a colocalized transcript through a lead-SNP or top-SNP shortcut. Unless the output identified a shared signal pair, the link between credible set and gene was reported as unresolved. Variant location within or near a lncRNA locus was not interpreted as evidence that the lncRNA molecule mediates the association.

A per-gene coding field was built from the reliable SuSiE credible-set members used in Figure 2D whose consequence was protein-altering (transcript ablation or amplification, stop gained or lost, start lost, frameshift, in-frame insertion or deletion, missense or protein-altering variant) or a canonical splice donor or acceptor variant. Each of the 147 qualifying coding variants was re-annotated with Ensembl VEP on GRCh37, the assembly of the GWAS positions [ref: VEP]. A gene was assigned when one of its protein-coding transcripts carried one of these terms. Gene IDs were mapped to GENCODE v49, by exact symbol when the GRCh37 ID was absent. Readthrough genes, genes without a symbol and genes whose GENCODE v49 gene type was not protein-coding were dropped. When a variant still hit more than one gene, we kept the gene whose term lay on the Ensembl canonical transcript. For each Catalog gene, the field records the maximum normalized SuSiE PIP (member PIP divided by the credible-set PIP sum) among its coding members, and the numbers of studies and sites. The unit is the gene as a Catalog row; studies are counted, not treated as replicates. No test was computed. The field misses non-canonical splice effects.

### Regulatory context of colocalized signals

For each colocalized gene–study pair, we selected one primary signal pair deterministically and intersected its normalized SNP.PP.H4 posterior probability with open chromatin from two liver snATAC sources across five shared lineages. Direct MASLD/liver-fat and liver-enzyme traits were kept separate. Counts use the gene as the unit, taking the best state over a gene's study pairs. Overlap below 5% of the posterior in both sources is grouped with absent overlap. The locus-matched expectation for a pair and lineage is the fraction of its mapped tested variants inside a peak of that lineage in either source. Observed and expected masses were averaged per gene and compared by 2,000 gene bootstrap resamples (seed 20260923). This expectation matches the window, its variant density and its peak coverage, but not LD or allele frequency. Partial or absent peak overlap is unresolved context, not a negative functional test.

For Figure 5E, a program member gene counts as open when one of its primary pairs places at least half of its posterior in open chromatin of the program's own lineage in at least one source (states replicated accessible or source-dependent). Hepatocyte, cholangiocyte, macrophage, fibroblast and T-cell programs map to hepatocyte, cholangiocyte, macrophage, stellate and T/NK chromatin. The count uses the same 767 colocalization pairs that Figure 5D summarizes per gene (22 direct-trait and 424 enzyme genes), including the three MVP European-ancestry studies.

### Single-cell atlas assembly and quality control

Seven public human liver single-cell or single-nucleus datasets were harmonized at donor and cell-type levels, comprising 104 descriptive biological donors, 275 integrated libraries and 1,232,318 cells [ref: source study for each dataset]. After two prespecified low-cell GSE189600 libraries were excluded, the analyzed atlas contains 102 donors, 273 libraries and 1,232,285 cells. GSE189600 deposits two sequencing runs for each of three healthy and three NASH replicate samples. The atlas ingests one run from four of them, and two libraries (one Healthy, one Steatohepatitis; 142 nuclei in the scored lineages) are analyzed. The deposit does not state whether the replicate labels are distinct people. Because ambient correction failed for this dataset, it is outside the primary stage-model universe and enters only the 62-donor sensitivity.

Per-library Cell Ranger filtered count matrices were integrated with ScaleSC on GPU [ref: ScaleSC]. Genes with fewer than three counts were removed, as were cells with fewer than 200 or more than 10,000 counts. We selected 4,000 highly variable genes, normalized and log-transformed the counts and computed 50 principal components. Harmony corrected the components with library as the batch variable (at most 30 iterations) [ref: Harmony]. Dataset labels were kept as metadata and were not the batch variable. A 30-nearest-neighbor graph on the 50 corrected components was used for Leiden clustering at resolution 0.5 and for UMAP [ref: Leiden] [ref: UMAP]. scVI was not used to integrate the atlas. Cell types were annotated with CellTypist and its Healthy_Human_Liver model with majority voting, applied to raw counts normalized to 10,000 per cell and log-transformed [ref: CellTypist]. The majority-voting label is the cell type used throughout (16 cell types), and the maximum prediction probability is kept as annotation confidence. We checked canonical marker expression per cell type; it did not change any label.

Program scoring used the integrated atlas with its original per-cell filters: a 200-count floor and a 10,000-count ceiling, no mitochondrial-fraction filter, no doublet removal and no ambient-RNA correction. A later harmonized QC (mitochondrial fraction < 20%, at least 200 genes and 500 counts, Scrublet score < 0.25) [ref: Scrublet] retains 1,030,882 cells. The share of scored-lineage cells that pass it is 49.4% in GSE185477, 79.3% in GSE244832, 89.1% in GSE174748, 97.0% in GSE202379 and 7.7% in GSE189600. Rescoring on passing cells keeps Ductular injury donor scores at Spearman 0.995 with the published scores.

Several datasets deposit multiple libraries per person that are not independent biological replicates: GSE244832 (117 sequencing runs from 18 donors), GSE202379 (67 from 46), GSE185477 (21 from 3) and GSE136103 (20 CD45-sorted libraries from 10 donors). The remaining datasets already had one library per donor. We therefore collapsed sequencing runs and CD45 sort-fraction libraries to biological donor before any per-donor statistic. Per-cell-type counts were summed and cells pooled to the donor before testing. Cell-type differential expression was tested by pseudobulk methods at the donor level, not by treating cells, runs or sort fractions as independent replicates. Where a disease-stage stratum fell below the prespecified minimum donor count, or was fully confounded with a single dataset, the contrast was reported as not estimable rather than significant. The atlas records the best cell-type effect, adjusted p value, number of significant cell types, hepatocyte effect and cell-type-specificity annotations.

### Single-cell diagnosis and stage labels

We corrected diagnosis labels against primary sources. The corrected labels cover the 102 analyzed donors: 24 Healthy, 14 Steatosis, 36 Steatohepatitis, 17 Cirrhosis and 11 Unresolved. The two GSE174748 NAFLD explants are explicitly cirrhotic in the Methods of Filliol et al. [ref: Filliol et al.; PMC9949942]. They receive ordinal code 3 and are excluded by the primary Healthy/Steatosis/Steatohepatitis estimand. Exact donor joins to Guilliams et al. Table S6 [ref: Guilliams et al., Table S6] give GSE192742 (Guilliams et al., 2022) 4 Healthy, 3 Steatosis, 1 Cirrhosis and 11 Unresolved donors, and all 19 remain excluded from stage analysis. Fibrosis without sufficient clinical information on etiology and activity was not recoded as steatohepatitis. GSE185477 Healthy labels are supported by the original description of the healthy-donor cohort [ref: GSE185477 healthy-donor cohort description; PMC8948611].

F-stage-specific single-cell analyses were restricted to GSE202379 (Gribben et al., 2024) samples with a documented SAF-derived F stage [ref: Gribben et al., 2024]. Documented F stage was recoverable for 58 of the 67 sequencing runs. Those 58 runs map to exactly 46 biological donors, with concordant stage labels across every multi-run donor. The accompanying scVI F-stage inference was validated at run level rather than donor level, so it does not meet the donor-level standard applied to every other per-donor statistic. We therefore excluded scVI-predicted and disease-stage-augmented F-stage fields from primary cross-cohort inference.

### Hotspot program discovery

The 117 Hotspot programs were discovered in the single-cell atlas [ref: Hotspot]. Their memberships, weights and names were fixed at discovery and were not changed or rediscovered by any analysis described here. Discovery kept a library's hepatocytes, cholangiocytes, fibroblasts or macrophages only when that library had at least 30 cells of the lineage, and kept T cells without a minimum. The program-stage models below score these same retained cells. Discovery ran separately in each lineage: hepatocytes (657,412 cells), fibroblasts (38,439), macrophages (19,595), cholangiocytes (32,402) and T cells (42,295). The fibroblast, macrophage, cholangiocyte and T-cell runs excluded GSE136103 and GSE192742. The hepatocyte run included every dataset.

Hotspot used its `danb` model on raw counts, with total counts per cell as the depth term. The latent space was a 30-component PCA of library-size-normalized, log-transformed counts for hepatocytes, fibroblasts, macrophages and cholangiocytes, and the atlas Harmony components for T cells. The neighbor graph was an unweighted 30-nearest-neighbor graph (50 for macrophages). Genes were restricted to protein-coding genes and lncRNAs. We removed genes on chrY and chrM, ribosomal (RPS, RPL, MRPS, MRPL), mitochondrial (MT-) and hemoglobin genes and a list of sex-chromosome genes; the four non-T-cell inputs also had a list of X-chromosome escape genes removed. Genes detected in fewer than 1% of the lineage's cells were dropped. Local correlations were computed for the 2,000 genes with the highest autocorrelation Z among genes with autocorrelation FDR < 0.05 and Z ≥ 7. Modules were called with `create_modules` (at least 20 genes, core genes only, FDR < 0.05), and each gene's weight is its mean local-correlation Z with the members of its module. This gave 30 hepatocyte, 29 fibroblast, 19 macrophage, 28 cholangiocyte and 11 T-cell programs (117 in total). A program's name records the cell-type subset it was discovered in.

Content lineage was assigned separately by scoring each program's fixed member genes against donor-level per-lineage pseudobulk. Scores were z-scored across lineages, so that the index reads specificity rather than abundance. The index was validated against donor-paired same-atlas contrasts before use. Content lineage is an added annotation and does not change program memberships, weights or names.

### Single-cell program stage models

The primary ordinal models use the 60 donors whose datasets could be ambient-corrected (GSE189600 excluded). Each program score was modeled as `score ~ stage + dataset` by generalized least squares with a separate residual variance per dataset (REML; `nlme::gls` with `varIdent`) [ref: nlme]. We used this model because the residual standard deviation differed between GSE244832 and GSE202379: the ratio is 3.6 for Ductular injury and 2.8 for Stromal ECM (variance ratios about 13 and 8), with a median of 1.6 and a range of 0.07 to 15.3 across the 117 programs. Donors in the two single-stage datasets (GSE174748 and GSE185477, Healthy only) carry no stage information and are not fitted. The p value is (1 + the number of 19,999 within-dataset stage permutations, seed 20260923, with |t| at least the observed value) / 20,000. BH correction is applied over all 117 programs, and selection also requires a stability median of at least 0.5. In simulation, the model-based 95% interval covered the true effect in 94.2% (true beta 0) and 95.5% (true beta 0.30) of draws. The 62-donor equal-variance OLS fit, with BH over 117 programs and an HC3 sensitivity, is a supplementary sensitivity (Supplementary Note 8).

Before the outcome was read, the family-level null (the first 2,000 of the 19,999 within-dataset stage permutations, BH over 117) gave at least one call in 4.5% of draws. We measured detection of each named program's earlier estimate (its beta in the 62-donor equal-variance fit) by adding that effect to the real scores with stage permuted within dataset (1,000 draws), and counting fits whose |t| exceeded the program's permutation critical value at 0.05/117. Detection was 0.999 for Ductular injury (earlier estimate 0.402) and 0.794 for Stromal ECM (0.352). A program below 0.80 is reported as not excluding its earlier effect. All checks passed, and no program was selected (0 of 117). The two named programs, Stromal ECM and Ductular injury, were selected from an earlier version of this stage model and kept fixed for the follow-up analyses below.

In the 62-donor equal-variance sensitivity, a 50-cell lineage gate removes the single Ductular injury call (q 0.0376 to 0.0765). Program gene coverage requires eight genes and 80% retained weight. Projected scores are not in native Hotspot units.

### Ambient RNA sensitivity

The atlas used for Hotspot discovery and primary program scoring did not receive ambient-RNA correction. In the aligned sensitivity, decontX was run against raw counts without changing the 117 prespecified program memberships, weights or labels [ref: decontX]. It succeeded in GSE136103, GSE174748, GSE185477, GSE202379 and GSE244832. GSE189600 and GSE192742 could not be corrected, because an ambient batch contained only one contamination cluster. Both are excluded from the aligned comparison and are kept only as explicit uncorrected passthrough datasets in a full-universe sensitivity. After the cirrhosis exclusion, the primary common universe contains at most 60 donors across four datasets. Corrected program scores were collapsed to donor and tested with model-based and HC3 standard errors. BH correction was applied separately over all 117 programs for raw, corrected and paired correction-induced effects. Integer-count transport was valid for 106 programs. The 11 prespecified T-cell programs were labeled untestable, because their native non-integer scoring data could not be reproduced. This sensitivity can qualify, but cannot establish, the source of ambient-sensitive transcripts. The corrected fit keeps all scores, the 117 memberships and weights, covariates and separate 117-program testing families; the eligibility rule removes only the two source-cirrhotic donors. Paired raw and corrected analyses use the same correction-available cells and donors, the original raw full-transcriptome exposures and the fixed raw gene calibration. Correction gaps and a separate full-raw analysis are reported explicitly.

### Content-lineage disease tests

Two targeted content-lineage tests used fixed projected weights in fibroblasts (Stromal ECM) and cholangiocytes (Ductular injury). Categorical steatohepatitis versus steatosis was primary. Steatohepatitis versus healthy was a separate secondary two-test BH family. Cohort-specific HC3 intervals and inverse-variance summaries required three donors per arm, a full-rank design and positive residual degrees of freedom. A 50-cell gate was primary, and 20 cells and annotation confidence ≥ 0.90 were sensitivities. Captured fractions among the six included lineages are not tissue proportions.

### Same-atlas localization of programs

We evaluated same-atlas localization of the two named programs with donor-paired, lineage-minus-hepatocyte score contrasts. Five lineages were prespecified for each program, giving a ten-contrast BH family with the same 50-cell gate. The focal comparison is the fixed lineage score paired within biological donor to that donor's hepatocyte score, with model-based and HC3 standard errors reported. This analysis describes cellular context within the same atlas and does not test transcript origin, secretion or cell autonomy.

### Independent single-nucleus datasets

We also scored the two named programs in an independent single-nucleus dataset (GSE256398), as hepatocyte pseudobulk from donors with at least 50 hepatocyte nuclei [ref: GSE256398 source study]. Two contrasts were prespecified: ordinal stage in non-cirrhotic donors (`score ~ stage + sex + age`; n = 13), and disease versus healthy (`score ~ disease + sex + age`; n = 17, with 6 healthy and 11 disease donors, 4 of them with MASH cirrhosis). Inference used HC3 standard errors with Holm correction over the four program-by-contrast tests.

GSE289173 (48 donors, single-nucleus RNA-seq) was opened to test whether the Ductular injury score is higher in later-stage donors [ref: GSE289173 source study]. The test could not be run. Only 14 of 21 program genes, carrying 0.64 of the fixed program weight, were detected in at least 1% of hepatocyte nuclei (prespecified minimum 0.80).

### Cell-cell communication

We inferred cell-cell communication with LIANA consensus analyses and follow-up methods suited to the modality, aggregated to biological donor [ref: LIANA]. Results are interpreted as inferred ligand-receptor coordination, not as direct measurements of signaling. The native input contained 117 units. The metadata join excluded 38 unmatched GSE192742 run aliases, leaving 79 units; this is distinct from the separate 125-unit lineage census. The source-qualified fit also excludes 10 GSE136103 donors under the stage-analysis exclusion, leaving 69 biological donors across four datasets. Among 55,402 ligand-receptor and cell-type groups, 349 models yield 1,047 stage coefficients, and 55,053 groups are untestable because of insufficient stage counts. BH correction was applied separately to 27 cell-type-pair-by-stage-term families, followed by Bonferroni correction across nine cell-type-pair families. Of the 13 display pairs, six are testable (all corrected q = 1) and seven are untestable.

For the supplementary program-linked communication screen, we crossed donor-level scores for the two named Hotspot programs with the 13-pair communication roster. For each donor and pair, the communication score was −log10 of the LIANA magnitude rank, lower-bounded at 10⁻⁴. Each of the 26 planned models used the communication score as the response and the standardized program score as the focal covariate, with disease stage, dataset, log10 sender cell count and log10 receiver cell count as covariates. Models with fewer than 30 matched biological donors were labeled untestable. Model-based and HC3 p values were BH-corrected separately over all 26 planned tests, and support required q < 0.05 in both families. Partial residuals were used only to display the adjusted supported coefficient. Leave-one-dataset-out directions were descriptive stability checks and did not define another significance family. Because program and communication scores come from the same atlas, we interpret their association as contextual covariation, not a signaling mechanism.

The screen has 22 estimable and four untestable tests within the 26-test family. Two Ductular injury associations pass both corrections, each with 58 donors: CDH1→PTPRM (beta = 0.4957; BH q = 3.628×10⁻⁶, HC3 q = 1.167×10⁻¹⁰) and ST6GAL1→EGFR (beta = -0.2987; q = 0.006914, HC3 q = 0.04880). Both keep their direction in four of four cohort omissions, with coefficients ranging from 0.003372 to 0.6983 and from -0.5662 to -0.06964, respectively. These are same-atlas RNA covariation results, not evidence of ligand secretion or receptor activation.

### Program scores in bulk tissue

For bulk transfer of tissue state, the prespecified L1-weighted program members were projected into the F1–F4-versus-F0 bulk fibrosis-stage models. Each displayed value is a weighted gene-level log2 fold change over the fixed members. Point area reports the fraction of available member weight that passes the complete 24,115-gene BH family within that stage contrast. These summaries are cross-sectional and descriptive; they are not a program-level significance test.

Figure 4E uses the 364 genetic-identity clusters of the three 774-unit cohorts that record both fibrosis and NAS (GSE130970 71, GSE135251 207, GSE162694 86). Per cohort, each of the 112 eligible fixed program scores was regressed by ordinary least squares on sex, centered fibrosis stage and centered NAS, with HC3 standard errors. The composition-adjusted model added three within-cohort principal components of centered log-ratio fractions and was refitted unchanged with each of two cell-type fraction sets (fragment-count and earlier-matrix BayesPrism fractions). Cohort estimates were combined by REML random-effects meta-analysis with Knapp-Hartung adjustment [ref: Knapp-Hartung], and BH correction covered 224 tests (112 programs by two axes). Adjusted calls had to be identical under both fraction sets. A negative control shuffled the earlier-set components within cohort 100 times (seed 20260923). A four-cohort sensitivity arm added 93 GSE174478 source records (457 records: GSE135251 207, GSE162694 86, GSE174478 93, GSE130970 71). These source-record counts are not verified independent-donor counts, and GSE174478 lies outside the 774-unit data.

Bulk transcription-factor activity was inferred separately for each F1–F4 versus F0 contrast with `decoupleR::run_wmean`, DoRothEA confidence A/B/C regulons, seed 42 and 1,000 permutations [ref: decoupleR] [ref: DoRothEA]. BH correction was applied across every tested TF within each contrast before HNF4A, RELA, SMAD3 and TCF7L2 were selected for display. Filled points require q < 0.05 in that complete contrast-specific family. The score represents inferred regulon activity, not TF expression or direct binding.

### Program-pair covariance in bulk liver and Visium tissue

*Discovery.* We computed bulk scores of the 113 programs that keep at least 80% of their fixed weight in bulk in four discovery cohorts (469 source records: GSE130970 76, GSE135251 214, GSE162694 86, GSE174478 93). Scores were residualized within cohort on sex, fibrosis stage and NAS bin (2 or less, 3 to 4, 5 or more). These source records were defined before the 774-unit genetic-identity roster and are not verified independent-donor counts. For each program pair, a partial correlation was computed after removing the leading residual component, with genes shared by the two programs left out of both scores. Of 6,195 pairs in the testing family, 2,426 had a positive partial correlation that passed BH correction for both the gene-disjoint and the partial statistic. Memberships, weights, the 2,426 pairs and four communities were fixed before the analyses below.

*Uncertainty of the fixed network.* We estimated uncertainty by resampling source participants within cohort, not edges, keeping all 113 programs and four communities. The summary for community C3 distinguishes within-content from cross-content pairs, and shared genes are removed subject to the 80% retained-weight rule. Histology terms are additive categorical fibrosis and NAS bins plus sex. Residual-score PCs were estimated per cohort and bootstrap draw, and the composition sensitivity added three centered log-ratio components. These pointwise intervals are conditional on the observed gene calibration and composition basis. Draws with a singular design, or with an observation whose hat value is 1, are reported as unusable rather than assigned a null effect.

*Two further bulk cohorts.* The unit was the participant, with one biopsy each. The primary sets were 91 of 94 GSE281797 participants (three removed because expressed sex contradicted the deposited sex) and the 87 of 93 GSE225740 participants with fibrosis recorded (sex inferred from expression) [ref: GSE281797 and GSE225740 source studies]. The secondary sets were participants staged F0-F1 (66 and 53). Counts were restricted to the discovery gene universe, TMM-normalized within each analysis set and scored with the same code under sex + fibrosis stage + NAS bin. GSE281797 kept 113 programs and 2,425 eligible coordinated pairs, and GSE225740 kept 110 programs and 2,320 pairs. The statistic was the fraction of coordinated pairs with a positive partial correlation after removing 1 residual component. It was also computed after removing 3 components, among the 1,448 pairs (1,368 in GSE225740) that were also positive at 3 components in discovery. The null maps each program to a random program of the same lineage (1,000 stored permutations), which keeps the lineage pair of every edge. The excess is the observed fraction minus the mean over 200 of these permutations. A participant bootstrap (1,000 draws within fibrosis stage) recomputed the excess end to end. A set's p is the largest of its permutation and bootstrap p values at both component counts, so it passes only if all four pass. Holm at family-wise 0.05 covers the two primary sets and, separately, the two F0-F1 sets, which were opened only if their cohort's primary set passed. The seed 20260923 was set before each random block.

Controls for this test were as follows. Our code reproduced the GSE193066 fraction (0.897) and comparator (0.709) to within 1e-9. The lineage null detected the GSE193066 signal (p = 0.000999). Gene permutation gave mean fractions of 0.485 and 0.483 (prespecified 0.44 to 0.54). The 869 and 820 pairs called informative-null in discovery were not above the null (p = 0.53 and 0.98). The label test rejected in 0 and 1 of 200 gene-permuted datasets (limit 16). The power check gave minimum detectable excesses at 1 component of 0.083 and 0.082, against the GSE193066 excess of 0.291. No numbered stop rule applied. A descriptive composition arm added three log-ratio components of InstaPrism proportions, log library size and top-100-gene count share to the design [ref: InstaPrism]. GSE281797 failed its acceptance rule (median hepatocyte fraction 0.44, accepted range 0.60 to 0.95) and was not read. Results are in Supplementary Note 10.

*GSE213621.* The unit was the genetic-identity cluster of the 774-unit roster. Coordinated pairs were rediscovered under a sex and three-level fibrosis design in 371 disease units of GSE130970 (65), GSE135251 (197) and GSE162694 (109), giving 2,328 pairs. These were to be tested in 258 GSE213621 disease units (87 F0-F1, 94 F2, 77 F3-F4) with four statistics: sign concordance, fixed-partition modularity against label permutation, the same modularity against rewiring, and ductular-portal community density. Holm correction covered the four statistics, and each had a unit-bootstrap condition (seed 20260923). The leave-GSE162694-out positive control passed (sign concordance 0.936 against 0.749). Two negative controls were run on permuted GSE213621 data before any primary statistic (100 permutations each). Each failed if its mean excess over the non-coordinated baseline exceeded 0.02, its 95th percentile reached 0.05, or its mean sign-matched excess exceeded 0.02. Gene permutation passed (−0.028, 0.010). Program-group permutation, which keeps covariance within each program's genes, failed (0.025, 0.081). Because a negative control failed, the prespecified stopping rule applied, and no GSE213621 outcome statistic was computed.

*Visium co-location.* The unit was the donor: 32 HMSMA Visium arrays, one per donor (9 female and 23 male), after three arrays from inferred controls were excluded [ref: HMSMA source study, HRA007511]. The source study lists 7 control, 10 MASL and 18 MASH livers without a key to the arrays, so 4 of the 32 retained donors are probably controls that cannot be identified. We scored per spot the 44 programs with at least 80% of their weight in genes detected in at least 1% of spots in every array (log-normalized counts, genes z-scored within array, fixed weights renormalized over detected genes). Scores were residualized on log library size and log genes detected, and Pearson correlations across spots were Fisher-z transformed within each array. The donor statistic is the mean z of coordinated minus non-coordinated cross-lineage pairs that share no genes (160 against 372; T1), and separately among pairs with positive bulk covariance (134 against 80; T2). The cohort statistic is the mean over donors. Each null permuted the coordinated label within lineage-pair strata (10,000 draws). T1 and T2 were one-sided, with Holm over the two at 0.05. Any reading beyond "not supported" also required that the per-donor T1 excess pass a sign-flip test and be positive in at least 22 of 32 donors, and that shuffling program identities within lineage (10,000 draws) not reproduce the T1 excess. A supported reading further required T2 to reject. Four secondary variants (zonation scores added as covariates, the first spatial component removed, gene-sharing pairs added with shared genes left out, program-specific terms removed) formed a Holm family of four that qualifies wording only. The seed was numpy SeedSequence(20260923) with 15 child streams. Controls: zonation scores reproduced their prespecified value (mean z −0.1705, negative in 31 of 32); a planted difference of 0.03 was detected in 0.998 of draws; random labels gave a T1 false-positive rate of 0.038; spot shuffling gave |D| of at most 0.0006; and portal-tract marker co-location was positive in 28 of 32. No stop rule applied.

### Single-nucleus ATAC-seq

The snATAC analyses used the two liver sources listed under Data availability (GSE244832 and GSE281367), which supplied lineage-resolved open chromatin [ref: GSE244832 and GSE281367 source studies]. In GSE244832, lineage labels were transferred from the donor-matched single-nucleus RNA atlas. ATAC gene-activity scores (normalized to 10,000 and log-transformed) were projected into each donor's RNA principal-component space (3,000 highly variable genes, 30 components), and each nucleus took the majority CellTypist label of its ten nearest RNA neighbors. In GSE281367, each nucleus received the label with the highest marker gene-activity score (summed activity of the marker genes divided by total activity) among ten marker sets. Labels were mapped to five lineages: hepatocyte, stellate (CellTypist fibroblasts in GSE244832), macrophage (including Kupffer cells in GSE281367), cholangiocyte and T/NK.

Peaks were called with MACS3 through SnapATAC2 in each cohort without condition labels (narrow peaks, shift −100, extension 200 bp, hg38 ENCODE blacklist v2 removed) [ref: SnapATAC2] [ref: MACS3]. A pooled lineage peak (q < 0.05) was kept when it overlapped a donor-level peak (q < 0.10; donor-lineage groups with at least 20 nuclei) in at least two donors. For each lineage, the peaks of the two cohorts were merged into fixed 500-bp intervals with SnapATAC2 `merge_peaks`, and intervals overlapping the blacklist were removed. Each count is one deduplicated fragment record overlapping the peak, summed in donor-by-lineage pseudobulks.

Differential accessibility compared MASH with normal donors that had at least 20 nuclei in the lineage. Peaks passed `edgeR::filterByExpr` with an intercept-only design in each cohort and were kept when testable in both cohorts. We then used TMM normalization, `limma::voomWithQualityWeights`, `lmFit` and robust `eBayes` with the design `~ group`, with BH correction within each lineage and cohort. A peak was supported when q < 0.05 in both cohorts with the same sign. Jointly testable peaks numbered 30,650 in hepatocytes and 9,264 in stellate cells, and no macrophage peak was jointly testable. With 5/9 and 6/6 normal/MASH donors per cohort, 38,066 of 39,914 jointly testable peaks are indeterminate and 1,844 are source-dependent (q < 0.05 in one cohort only; 1,505 share sign with the other cohort).

### Spatial transcriptomic analyses

We integrated spatial organization, zonation and cross-cohort annotations while keeping each source's measured-gene background, biological unit and provenance. None of the spatial analyses supplies a valid MASLD-versus-healthy disease-effect estimand. GSE192741 and Vu et al. [ref: Vu et al.] test residual spatial autocorrelation against matched genes. Yakubovsky et al. is a non-MASLD tissue-context analysis [ref: Yakubovsky et al., 2026]. GSE287826 was not analyzed, because no authoritative donor key was public. Spatial results therefore localize fixed genes and prespecified programs in tissue; we do not describe them as independent disease-associated change.

Visium inputs were the source-processed Space Ranger outputs, without realignment: five human GSE192741 sections (mouse sections were not used) and ten Vu et al. CytAssist FFPE arrays. Per Guilliams et al. Table S6, the four GSE192741 donors (H35–H38) were sampled at gastric bypass (H35, H37), cholecystectomy (H38) or resection of colorectal-cancer liver metastases (H36, who reported six standard alcohol units daily), with recorded steatosis of 35–40%, none, 70% and 10% and fibrosis of none, mild, clear and mild. H37 and H38 also contributed single-nucleus libraries to the atlas. Spots were kept with at least 200 detected genes, at least 500 counts and at most 20% mitochondrial counts, and genes detected in at least ten spots were kept. Cell-type abundance per spot was estimated with cell2location [ref: cell2location]. The reference was a negative binomial regression model trained on raw counts of the CellTypist-annotated single-cell atlas, with the 16 cell types as labels and dataset as batch. It used a different build of the atlas (1,170,457 cells from six datasets), subsampled to 150,000 cells with an equal cap per cell type (seed 42), genes detected in at least ten cells, 250 epochs, batch size 2,500, learning rate 0.002 and 1,000 posterior samples. The spatial model used the genes shared with the reference and detected in at least five spots (10,738 for GSE192741 and 10,056 for Vu et al.), section or array as batch, 15 expected cells per spot and detection alpha 20. It trained for 30,000 epochs on the 6,546 GSE192741 spots and 15,000 on Vu et al., with 1,000 posterior samples. Downstream analyses use the 5% posterior quantile of abundance (q05).

Each Visium spot received a program score: the z-scored expression of its measured fixed member genes, weighted by their original L1 weights rescaled to sum to one. The prespecified confirmatory design residualizes the score for the cell2location hepatocyte abundance (the discovery lineage of both programs), log1p library size and log1p detected genes. A model residualized for all 16 cell2location q05 abundance factors and the same two covariates was added afterwards as a supplementary falsification test of composition dependence; Figure 5F plots it. Residual Moran's I on six-neighbor spot graphs built within connected tissue islands [ref: Moran's I] was compared with 9,999 random gene sets scored and residualized identically (seed 42 for GSE192741, 1042 for Vu). Random sets were matched on expression, detection rate, biotype, mitochondrial and ribosomal status and correlation with the hepatocyte abundance. Expression, detection and correlation bins were relaxed by up to two bins when no match was found. The reported excess is the observed statistic minus the null mean. p = (k + 1)/(9,999 + 1), where k counts draws at or above the observed value, so the smallest attainable p is 1e-4. BH correction covers the two prespecified programs within each source; no prespecified test compares the two programs. The null resamples genes and holds donors and arrays fixed. Donor-first estimates are used where a donor key exists (five GSE192741 sections collapse to four donors). An exact sign-flip test of per-donor excess has a smallest attainable p of 0.0625 with four donors and is reported descriptively. Array-level results without a donor key are reported as source-dependent and never as donor-level inference. Positive spatial autocorrelation without matched calibration does not establish program-specific tissue architecture.

### Assay observability of programs

Complete per-assay coverage is in [supplementary observability table; number to be assigned]. CosMx adequately observes 12 of 117 programs, because its panel is targeted. GeoMx was excluded from inference because no gene axis was loaded and its direction and normalization were unresolved. GSE287826 (no gene axis loaded; no AOI-to-donor key) tests no complete program in the current processing, so its observability is not yet determined. The four Visium/ST datasets each observe 114 of 117 programs. The Gene Catalog defines no negative-call rule, so a tested-negative state is not assessed. These are limits on what the assembled assays and data can resolve, and they are the reason that counts of supporting modalities are not used anywhere in this work.

### Proteomic analyses

The proteomic endpoint uses the PXD051911 liver DIA-MS cohort (58 participants) [ref: PXD051911 source study]. PXD052937 plasma was excluded because its per-sample disease labels could not be authenticated. Protein significance, effect direction, measured-gene coverage, cohort overlap and selection conditioning were recorded explicitly. We describe protein results as tissue context or corroboration in the native assay, not as functional validation. Analyses started from the liver protein-group intensity table deposited by the source study (7,096 protein groups in 58 runs, on a linear scale, with missing values left missing). The source study performed spectral search, protein quantification and any normalization before deposit, and we did not repeat these steps.

Differential abundance was fitted with `~ group + acquisition_batch + age + BMI + sex`, and BH correction covered the complete 6,071-protein family. The 25-protein process view in Figure 5C was selected using protein significance and RNA-protein direction concordance in this same cohort. It is therefore descriptive and is not an independent validation of the bulk or single-cell programs.

For Figure 5E, single-gene protein groups with at most 50% missing values were log2-transformed, quantile-normalized, collapsed to one value per symbol by the median and z-scored per protein across the 58 samples. Each of the 117 fixed programs was scored per participant as the L1-weighted mean of its measured member proteins, with members counted once per primary-assembly gene (the rule used for the gene column of Figure 5E). A program was scorable with at least eight measured members and at least 20% of its L1 weight retained (95 of 117). Histologic burden was the first principal component of steatosis, ballooning and lobular inflammation grades and Kleiner fibrosis stage in the 46 MASLD participants with complete histology. NAS was excluded because it sums three of these. The association was a rank-residual partial Spearman correlation adjusted for acquisition batch, age, BMI and sex. Inference used a batch-stratified Freedman-Lane permutation p (9,999 permutations), a participant bootstrap 95% CI resampled within batch (4,999 draws) and BH correction over the 95 scorable programs. A program was robust when q < 0.05 and the sign agreed for the equal-weight score, the score without its top-weighted gene and both acquisition batches. Each program's seeds were a base seed (20260923 for permutation, 30260923 for bootstrap) plus a fixed program index. The index was the program's position (1 to 94) in the program list among the 94 programs scorable under the earlier gene-symbol membership, or 1032 for Portal fibroblast (NR1H4), which is scorable only under the primary-assembly membership.

For the GNMT–MAT1A–CYP2C19 analysis, protein abundances in 46 MASLD participants were rank transformed and residualized for acquisition batch, age, BMI and sex. Pairwise Spearman correlations were computed between the residuals, with BH correction over the three prespecified gene pairs. Sensitivities added fibrosis stage, NAS or the four component histology scores to the residualization model, and correction was reported both within each three-pair family and across all 12 pair-by-adjustment tests. The protein covariance was reconstructed from PXD051911 participant-level inputs with the same 46 complete cases across the baseline, fibrosis, NAS and component adjustments. Ranked protein residual correlations used the degrees of freedom of the nuisance model and 2,000 paired resamples. A broader metabolic-abundance sensitivity fixed the union of three Hallmark metabolic gene sets before comparative results were inspected. All protein-trio covariance is conditional on how the trio was selected.

### External cross-reference with independent liver studies

We used external liver protein-QTL and single-cell multiomic studies only to cross-reference our evidence classes against source-reported prioritized genes, loci and transcription factors. Neither study contributed to class definition. These studies were treated as cited comparisons, not reanalyzed data. We compared their reported calls with our per-gene values, reporting direct-MASLD and enzyme-trait anchors separately. Loci where layers colocalized but effect directions differed were reported as shared loci, not as sign-concordant validations.

For Elison et al. [ref: Elison et al., preprint], the genome-wide supplementary tables of cell-type quantitative-trait, gene-regulatory-network and chromatin-loop statistics were embargoed at the time of analysis. We therefore cross-referenced the per-locus and per-transcription-factor calls stated in the preprint main text and figure legends: the featured colocalized loci, the cell-type-annotated risk loci and the stage-specific hepatocyte gene-regulatory-network transcription factors. No primary data from either study were reanalyzed or redistributed, and both are cited as preprints that have not been peer reviewed. Each comparison is recorded as external annotation of the evidence-class table and can be extended to the full supplementary statistics when they are released.

### MASLD Gene Catalog construction

Gene Catalog entries were built only from checksummed output files of the analyses described above. Each gene or program keeps a stable identifier, evidence class, testing universe, assay-specific call, source dataset, relationship to the source publication and the graph of parent files. The evidence rows carry the native estimate and, for RNA, the BH q, but no standard error, interval, p value or unit count. Membership in a program is not expanded into gene-level protein, spatial, chromatin or functional evidence.

For a given assay, a gene, program, peak or variant receives one of seven evidence states: supported, discordant, tested negative, indeterminate, untestable, not applicable or source-dependent. A blank or non-significant comparison may reflect absent assay coverage, insufficient expression, weak molecular-QTL power, phenotype mismatch, an unavailable donor key, reuse of a source or a mechanism outside the resolution of the assay. A tested-negative state is reserved for an adequately observed, applicable and properly joined assay that fails a prespecified adequate-negative or equivalence rule. Non-significance, inadequate coverage, a failed source check or an assay that cannot observe the proposed mechanism remains indeterminate or untestable. The Gene Catalog defines no such rule for either domain, so no gene can be tested negative. A gene outside the eligible participant or measurement universe of a contrast is untestable for that contrast, not a negative result. The states describe the role of the evidence. They do not assign cause, give a probability that a target is valid, or rank genes on the assumption that agreement across modalities is best.

### Follow-up experiment routes

Next-experiment routes are deterministic functions of the evidence, molecular object and testability fields. A noncoding regulatory element routes to allele-aware reporter assays or locus editing in the lineage where the element is accessible. A disease-associated lncRNA without a genetic anchor routes to RNA depletion with readouts of neighboring genes and phenotype. A colocalized lncRNA or ambiguous lncRNA locus routes to paired RNA depletion and locus-level CRISPRi or allele perturbation, so that effects of the RNA product, of transcription and of the DNA element can be distinguished. Protein-coding genes keep the activation, inhibition, protein-level, lineage-restricted and context-specific routes. The routing does not select new scientific results and does not use downstream outcomes. The signed Catalog fields are provisional, because the bulk, stage, composition and colocalization results they combine were not all computed on the same version of the data.

### Gene Catalog annotation fields

Each field below was specified, with its controls and stop rules, before its outcome was computed, and its states were assigned in code by the prespecified reading rule. The seed was 20260923. Class-dependent counts use the Catalog evidence classes, the colocalization tables and the 774-unit bulk fit. No field enters an evidence class or a next-experiment route. The coding-variant, caQTL and postmortem-fat fields, and the composition sensitivity that added no field, are described in their own subsections above.

*Healthy-liver zonation.* For each Catalog gene, we computed mRNA position along the central-portal axis from the Yakubovsky et al. (2026) Visium deposit in five non-steatotic live donors (M4 to M8; unit the donor) [ref: Yakubovsky et al., 2026]. We used spots in the source's zones 1 (central) to 8 (portal) and divided expression by spot total. Position was the center of mass of the zone-mean profile minus 4.5, taken as the median over the five donors. A gene was called central or portal only if it was detected in at least 5% of spots in every donor, fell on the same side in all five donors and passed a detection-decile threshold set from a zone-label permutation null (200 permutations). Genes below the detection rule are described as not measurable in Visium, not as unzonated. Hepatocyte protein zonation was copied from Weiss et al. (2026) (14 healthy individuals, Supplementary Table 3) and was not refitted [ref: Weiss et al., 2026].

The one prespecified test was the Spearman correlation between mRNA position and the Weiss coefficient over the 1,337 Catalog genes measured in both layers, excluding the landmark genes that define the zone axis. It used a two-sided gene-permutation p (10,000 permutations; family of one, alpha 0.05) and a 2,000-draw gene-bootstrap interval. A combined layer flag required p below 0.05, rho of at least 0.30 and a positive per-donor correlation in each of the five donors. Because genes covary, the permutation p is too small, so the donor condition and the rho floor govern the reading. The seed was numpy SeedSequence(20260923) with four streams. Controls: zone orientation was correct in 5 of 5 donors; 19 of 22 robust landmark genes had the same sign in both layers (rule at least 15); and 1.5% of zone-label permutations gave |rho| of 0.08 or more (rule at most 10%). No stop rule applied. Across the 1,337 genes measured in both layers, mRNA and protein positions agreed in direction (Spearman rho 0.39, 95% interval 0.34 to 0.44; positive in each of five donors). Weiss et al. reported a moderate RNA-protein zonation correlation against single-nucleus RNA (their Supplementary Fig. 9). This annotation describes healthy liver only and gives no information about zonation in MASLD.

*Assay observability.* Protein detection was taken from the PXD051911 proteomes: liver, 58 participants; plasma, 184 samples (143 baseline and 41 follow-up) from 143 participants. A gene was detected if its protein was observed in at least one sample. The field also records detection in at least half of the participants, using baseline samples only for plasma. CosMx targeting was read from the panel header of the GSE312698 CosMx data (Boesch et al., 2026) [ref: Boesch et al., 2026]: 942 Catalog genes have a single-gene probe, and 40 are covered by 19 probes shared between genes. Single-cell detection used donor pseudobulk from 104 donors (275 runs) in nine lineages. A gene was detected in a donor-lineage with at least 5 counts and at least 10 CPM, and in a lineage when detected in at least half of the eligible donors. A donor was eligible with at least 20 cells and 10,000 counts in that lineage (25 to 95 eligible donors per lineage).

Controls were known liver and plasma proteins and lineage markers (liver 7 of 7, plasma 8 of 8, single cell 10 of 11 listed markers); olfactory-receptor genes as negatives (1 of 387 detected in the liver proteome, 0 of 387 in plasma and on the CosMx panel, 0 of 292 measurable in single cell); and a rise in detection with bulk RNA abundance (highest minus lowest decile at least 0.30; liver 0.815, single cell 0.727, plasma 0.123). A field failing the known-answer, negative or coverage check would not be released. A field failing the abundance check is released with a note, which applied to plasma. There is no confirmatory test. Absence of detection is recorded as untestable, not negative.

*Liver protein coverage by class.* A gene was counted as quantified when the PXD051911 liver differential table (6,071 proteins; 58 participants, 12 control and 46 MASLD) contained a row for it after symbol updating with HGNChelper (90 symbols remapped) [ref: HGNChelper]. Coverage is a count by Catalog class, and no test was performed. Controls: 14 of 16 known liver proteins were quantified (PNPLA3 and TM6SF2 absent), 3 of 10 tissue-restricted non-liver proteins were quantified, and 89 of 89 remapped symbols were confirmed by an independent alias table. Detection in any sample (previous field) is a broader criterion than a row in the differential table.

*Follow-up annotation.* For every Catalog gene, we added gnomAD v4.1 LOEUF and pLI for the Ensembl MANE Select transcript (deciles over 17,325 genes) [ref: gnomAD v4.1], the DepMap Public 24Q4 common-essential call and median Chronos gene effect [ref: DepMap], and the Human Protein Atlas secretome location and blood concentrations [ref: Human Protein Atlas]. These fields describe how a gene might behave in a follow-up experiment and do not enter any evidence class or route. Controls used core-essential genes (CEGv2) and nonessential genes (NEGv1). The lowest LOEUF decile held 0.285 of core-essential genes against 0.024 of nonessential genes, and 0.914 of core-essential against 0.006 of nonessential genes were DepMap common-essential. All 10 known blood-secreted proteins tested were annotated as secreted to blood. No test was performed.

### Predictive-model evaluation

*Single-cell encoders.* The single-cell encoder comparison maps five broad cell classes, not disease stages. TranscriptFormer TF-Sapiens exposure is tied to the exact evaluated checkpoint (SHA-256 eff027e7393f92a32cadcbb2e13a917e24b62b71341571bf40f38ff6ec6d3021) through all 30 prediction-shard records and the 50,000-cell embedding record [ref: TranscriptFormer]. Its complete inventory of 644 human training datasets shows five Resource studies absent from training and two seen in training (GSE185477 and GSE192742). Pooled scores and displayed gains keep all seven studies and cannot be read as evidence from unseen studies only. Geneformer V2 has no qualified unseen study among those seven [ref: Geneformer]. The TranscriptFormer comparison restricted to unseen studies uses the original equal-study mean over five studies and 80 biological donor identifiers. Its 10,000 saved paired draws resample donors within study (seed 20260824) and preserve evaluable class coverage; cells are not resampled as independent units. Intervals are conditional on the fixed models, five studies, class coverage and observed cells, and do not estimate generalization to a population of studies. No models or scores were refitted for this comparison.

*Severity model populations.* Model-specific training and evaluation populations are reported separately from the population of the Resource differential-expression analysis. The released severity model card records 521 training sample records and no complete donor key. We do not describe these as 521 verified participants or as the bulk differential-expression roster. The released scorer applies its stored expression transformation and weights without refitting to evaluation outcomes. Its GSE268273 evaluation uses 109 participant records. Figure 7D evaluates GSE213621 in 263 genetic identities (one record each, after removal of 11 identities whose two libraries carry different stage labels) and GSE276114 as 81 MASLD and 82 chronic viral hepatitis participants. The published record-level contrasts use 299 disease-only and 367 control-inclusive GSE213621 records and 177 GSE276114 records. These are two additional cohorts and three contrasts, not three independent cohorts. Source-declared participant keys and joins of counts to phenotypes were checked; genetic identity was not verified.

*Exposure of evaluation cohorts.* Whether a cohort is external depends on the model. GSE213621 also contributes to Resource analyses, and all 109 GSE268273 records occur in BulkFormer pretraining [ref: BulkFormer]. The latter limits an unseen-sample claim for BulkFormer, not for the separately fitted severity scorer. Cohorts used to inspect refusal thresholds or other acceptance criteria are exposed to development for those criteria, and GSE268273 is exposed to development for the choice of expression transform. Four arms were eligible: each has a coherent expression transform and can score one sample at a time (out-of-cohort Spearman 0.6506 to 0.6631 on these 109 records). Among them, the released scorer is the one with the highest Spearman, a rule from the original prespecification. We amended the protocol to make coherence and reliability preconditions for eligibility after results for the first three transform families had been seen and before the quantile-normalization family, which the released scorer uses, had been run. Input coverage, score refusals and diagnostics are reported separately, and a diagnostic threshold not shipped with the scorer is not described as an enforced refusal rule.

*Comparison with in-cohort models.* We compared severity transfer with the recorded in-cohort ridge on the same evaluation population. In the three additional contrasts, the native comparator used full-cohort expression scaling before outer-fold evaluation, although feature selection and PCA were fitted on outer-training records. Those preprocessing steps were not refitted within each inner tuning split, so this comparator is transductive rather than a strictly inductive nested-validation estimate. The 10,000-draw stage-stratified paired participant bootstrap used seed 213621. Its intervals are conditional on the fitted predictions and do not include uncertainty from retraining. Participant-level predictions of the native comparator were not saved, so these paired intervals are source-reported rather than independently recomputed. An interval that includes zero does not establish equivalence.

*Agreement metrics.* Spearman correlation measures ordering, and its ratio to a tied-label ceiling is labeled separately from the raw correlation. Neither establishes agreement with histological stage or calibration on an absolute stage scale. The repeated-biopsy ICC comparison uses 58 source-declared GSE193066 pairs whose identities are not verified. It is a source-reported diagnostic, not primary evidence of within-person reliability or a calibration of label error. Possible biological change and shared technical effects further limit it, and a larger ICC alone does not establish greater accuracy. The model card reports that only three of these pairs have both biopsies within the released scorer's acceptance bounds. The 58-pair ICC therefore describes scoring outside the scorer's stated input limits, not performance on accepted inputs, and it was not recalculated here. Its stage SD of 0.6433 and the separately reported PXD051911 value of 0.4655 are inputs to sensitivity scenarios. They are not an empirically qualified 0.4655–0.6433 label-error range or a calibrated severity ceiling.

*RNA to H3K27ac.* For RNA–H3K27ac prediction, we checked the 99 participant keys and five fixed folds. The source rank-4 model uses training-fold preprocessing and a baseline equal to the training-target mean. Its recorded squared-error skill and uncertainty were not recomputed from participant losses. The forward and reverse skills are measured in different target spaces, so they cannot be compared to rank the models.

*Reporter-assay rows.* The Figure 7B reporter-assay rows score 1,033 GSE281364 elements in the HepG2 control and palmitate-oleate contexts only (tanh of the mean Fisher-z Spearman over the two contexts; 10,000 resamples of 239 long-range blocks) [ref: GSE281364 source study]. LX-2 control and TGF-beta outcomes are in the same outcome table but were not scored. The comparator is a ridge on allele identity. deltaSVM, gkm-SVM, sequence CNN or transformer baselines and MPRALegNet were not run. The Bonferroni family is the seven displayed reporter-assay rows, and the corrected lower bound is the 0.025/7 quantile of the same paired bootstrap draws. HyenaDNA's corrected lower bound stays above zero (+0.006) when the family is widened to all eleven fitted pretrained heads [ref: HyenaDNA].

*Scope.* Of the transfer metrics, four were reproduced from saved outputs; the others are source-reported comparisons or unresolved checks [model-evaluation summary table; number to be assigned]. No new model was trained for this assessment. Linear regression can use log-expression predictors; the transformation changes the relationship modeled and does not invalidate regression. Results for one scVI task do not imply that linear models are generally superior. A separate exploratory development (a rank-12 RNA-to-H3K27ac model with local correction, its acquisition and input-robustness checks, and an AlphaGenome allele analysis against Currin caQTL effects) is not part of Figure 7, and no result from it is reported in this paper.

## Quantification and statistical analysis

### Statistics and reproducibility

Biological donors or participants are the unit of inference in every test, except where a test is explicitly at the level of genes, loci or credible sets. Those units are comparison units, not replicates (see Units of inference). Each test reports its analysis universe, numerator, denominator, effect size and uncertainty where available. BH adjustment was applied within prespecified analysis families, and permutation p values used the plus-one correction. Null results and failed acceptance criteria were kept and used to limit the wording of claims. Nominal p values are never described as significant.

Primary definitions are distinguished from sensitivity analyses. Primary transcriptomic evidence is the pooled limma-voom quality-weighted result at BH adjusted p < 0.05 and |log2FC| > 0.50, with the effect-size-aware interval-null FDR gate kept as a comparator. The genetic analysis has two coverage layers across the 35 prespecified tier-1/2 GWAS strata: ABF colocalization for every testable gene–GWAS pair, and multi-signal SuSiE colocalization only where its prerequisites are met. The evidence classes use the SuSiE call, and the ABF map and the SuSiE-or-ABF union are reported with their own denominators.

The testing families stated with each analysis are collected below.

| Analysis | Family and correction |
|---|---|
| Disease-versus-control bulk differential expression | BH over 24,115 genes; interval-null comparator BH over the same 24,115 genes |
| GSE126848 within-cohort refit | BH over 15,592 genes (a separate family) |
| Fibrosis-stage and NAS contrasts; F1 and F4 cohort interaction | BH over 24,115 genes per contrast; 2-df moderated F test with BH over 24,115 genes |
| Library properties of stage contrasts | F1 versus F0 a family of one; other ten contrasts one BH family |
| Grouped-stage comparison in GSE213621 | Holm over two contrasts, family-wise alpha 0.05 |
| GTEx postmortem fat | Primary test a family of one (two-sided); two qualifiers a Holm family of two |
| Adolescent cohort (GSE66676) | One prespecified test, no correction |
| Genetic ancestry sensitivity | Holm over F1 = {T1, T2}, family-wise alpha 0.05 |
| Continuum gene-level models | BH over 24,115 genes; signature genes BH over 143 genes |
| Continuum pathway and program models | BH within each collection per score and model (50 Hallmark, 658 KEGG, 1,787 Reactome, 7,583 GO BP, 1,855 GO MF, 1,042 GO CC, 268 TF, 117 Hotspot, ten NMF) |
| Within-stage program models | BH over 117 programs per cohort, axis and model, and per meta-analysis axis and model |
| Continuum score sensitivities | 117 programs and 43 systems as separate BH families; ordinary and HC3 separate |
| caQTL comparison | One prespecified test, family of one |
| Figure 1B overlap | Exact stratified hypergeometric test, prespecified |
| Kumak et al. list comparison | Holm over two lists, family-wise alpha 0.05 |
| TWAS-to-disease comparison | Family of one, two-sided alpha 0.05 |
| Coding-mass differences | Bonferroni over three members (98.333% intervals); noncoding contexts Bonferroni over five (99% intervals) |
| Single-cell program stage models | BH over 117 programs; ambient sensitivity BH over 117 programs separately for raw, corrected and paired effects |
| Content-lineage tests | Secondary contrast a two-test BH family |
| Same-atlas localization | BH over ten contrasts |
| GSE256398 | Holm over four program-by-contrast tests |
| Cell-cell communication | BH within 27 cell-type-pair-by-stage-term families, then Bonferroni across nine cell-type-pair families |
| Program-linked communication screen | BH over 26 planned tests, model-based and HC3 separately |
| Figure 4E program models, unadjusted and composition-adjusted | BH over 224 tests (112 programs by two axes) |
| Bulk TF activity | BH over all tested TFs within each contrast |
| Program-pair covariance | BH over 6,195 pairs in discovery; Holm over two primary sets and separately two F0-F1 sets; GSE213621 Holm over four statistics |
| Visium co-location | Holm over T1 and T2; secondary variants a Holm family of four |
| Spatial matched nulls | BH over two programs within each source |
| Liver proteome | BH over 6,071 proteins; Figure 5E BH over 95 scorable programs |
| GNMT–MAT1A–CYP2C19 protein covariance | BH over three pairs; and across 12 pair-by-adjustment tests |
| Healthy-liver zonation | One prespecified test, family of one, alpha 0.05 |
| Reporter-assay rows (Figure 7B) | Bonferroni over seven rows |

Spike-in calibration, the descriptive molecular-systems summaries, the observability and coverage fields and the follow-up annotations have no hypothesis test and no family.

Seeds: 20260923 for most permutation, bootstrap and simulation analyses, as stated with each; 20260929 for the genetic-identity sensitivity; 20260907 for the class-comparison sign-flip diagnostic; 20260906 for the continuum extension; 20260923 + r for spike-in replicate r and 20261923 for its bootstrap; 30260923 plus program position for the proteome bootstrap; 42 for TF activity and for GSE192741 spatial nulls; 1042 for Vu spatial nulls; 20260824 for the TranscriptFormer donor bootstrap; 213621 for the severity-transfer bootstrap; 20260817 for the pooled Figure 2D interval; 20261017 for the Figure S2S subclass and Figure 4G and 4H system bootstraps; 9873 for the encoder donor-cluster bootstrap in Figure 7C and 7F (2,000 draws); and 20260831 for the LD-block bootstrap of sequence-model gains. Every randomized analysis records its seed and complete permutation family.

# Supplementary Notes

## Supplementary Note 1. Scope of the genetic and disease-state evidence layers

Liver-trait genetic evidence in this Resource takes three forms: multi-signal liver eQTL colocalization, single-signal liver eQTL colocalization, and fine-mapped protein-altering or canonical-splice credible-set members. The Gene Catalog states which form applies to each gene. None of them alone identifies the gene that mediates a trait association. Of the 436 primary multi-signal colocalized genes, 423 take their strongest colocalization evidence from a liver-enzyme trait and 13 from a direct disease or liver-fat trait.

Disease-state liver expression has two parts. The first is the pooled, cohort-adjusted disease-versus-control calls from 774 genetic-identity units in five cohorts (138 control, 636 disease). The second is the cross-sectional fibrosis-stage contrasts. The disease-versus-control fit tested 24,115 genes. The units are genetic clusters, not verified donors. Stage contrasts compare different people, not change within a person.

We also describe how fine-mapped posterior probability at liver-trait loci divides among protein-altering, splice and noncoding variants. In that analysis, the contrast between direct traits and liver enzymes is descriptive. It rests on 19 direct-trait regions and is pooled across ancestries. It is confounded with GWAS ancestry and fine-mapping coverage, and no ancestry stratum supports a trait-class call. It does not resolve the responsible variant at any locus.

## Supplementary Note 2. GNMT colocalization and linkage, and stage association of the two focal programs

GNMT has a single-signal colocalization with direct NAFLD (ABF PP.H4 0.706) but no multi-signal SuSiE call (best PP.H4 0.128). The disease side is sub-threshold (anchor GWAS p = 2.6 × 10⁻⁴; smallest p within 1 Mb 1.4 × 10⁻⁴), and PP.H4 depends on the prior p12: 0.33, 0.71 and 0.83 at p12 = 1 × 10⁻⁶, 5 × 10⁻⁶ and 1 × 10⁻⁵ in the earlier run of this pair. Its Catalog genetics call is therefore indeterminate, and its local regulatory haplotype is unresolved.

The rs2296804 anchor has SNP.PP.H4 0.259. rs2296805 has UK Biobank European r2 = 0.998 with rs2296804, and rs9471976 and rs11752813 have r2 = 0.818 and 0.580 with rs2296805.

For Stromal ECM (IGFBP7), the bulk fibrosis association does not reach BH q < 0.05 in the 364 units from the three cohorts of the 774-unit set that record both fibrosis and NAS (GSE130970, GSE135251 and GSE162694) (q = 0.214; q = 0.421 after composition adjustment). It reaches q = 0.045 only when 93 GSE174478 source records are added and the Knapp-Hartung variance is not floored at the REML value. In that analysis, q is 0.230 after composition adjustment. In the single-cell stage models, once dataset-specific variances are modeled, neither focal program passes within-dataset permutation BH over all 117 programs for a cross-sectional stage association (smallest q of any program 0.323).

## Supplementary Note 3. Gene-level comparison of colocalized and disease-state genes

### Overlap of the two gene sets

Multi-signal colocalization evaluated 3,769 genes of the native RNA–genetics interface. Of these, 426 were colocalized and 220 were canonical disease-state calls. In all, 26 genes were in both, against 21.21 expected within expression deciles (Mantel-Haenszel odds ratio 1.295, 95% CI 0.835 to 2.007; exact stratified *p* = 0.2755; one prespecified test). Rounded, these values are 21.2 expected and an odds ratio of 1.29 (95% CI 0.84 to 2.01). The deciles did not fully control expression, so we read no direction from this result.

The overlap is not informative in itself. The interval includes 1. The unstratified test had less than 80% power to detect odds ratios between about 0.46 and 1.77. Two sets of these sizes are mostly disjoint whether or not the biology differs. Among the 14,704 genes that both analyses list, 935 of the 1,129 disease-state calls without colocalization had no multi-signal colocalization result. These genes are untested by multi-signal colocalization, not genetic negatives. The annotation-adjudicated Gene Catalog has a separate 14,682-gene interface, and we did not substitute it for the universe of the native RNA comparison.

### Expression-matched comparison of hepatocyte specificity

Matched for expression, colocalized genes are not more hepatocyte-specific than disease-state genes, so the prespecified directional prediction failed. The comparison contains 281 expression-matched gene pairs (expression SMD 0.00371). In a post hoc comparison, the mean genetic-minus-state difference went the other way (−0.265515 log2 units; 126 positive and 155 negative pairs). The relabel-and-rematch outputs (*p* = 0.00079992; BH *q* = 0.00159984 across two outcomes) are diagnostic only. That null destroys the observed expression imbalance and matching structure, so it is anti-conservative, and family correction cannot repair it.

The separate post hoc paired sign-flip diagnostic gives nominal *p* = 0.01482985 (100,000 draws, seed 20260907). It assumes exchangeability within pairs. Flipping signs jointly within the 216 loci that hold the 281 colocalized genes gives nominal *p* = 0.01606 (100,000 draws, seed 20260907; 0.01625 with seed 20260923). A sign test gives *p* = 0.094677 (126 of 281 pairs positive), and 99 of 216 locus means are positive (*p* = 0.247). All of these tests are post hoc, and none is calibrated. The result is a descriptive difference conditional on selection. It is not confirmatory class separation or independent biological replication.

Of the 1,129 state-only genes, 1,124 have a single-cell hepatocyte-specificity value and form the comparison pool. Of these 1,124, 931 had no multi-signal colocalization result, so most of the pool consists of genes that colocalization could not test. The 193 tested genes have mean AveExpr 1.3405, against 0.4164 for the whole pool. In a second post hoc comparison restricted to those 193 genes, the difference kept its direction. In 129 expression-matched pairs (expression SMD −0.0004), the mean genetic-minus-state difference was −0.342 log2 units (median −0.183, paired d −0.217; nominal sign-flip *p* = 0.0154, 100,000 draws, seed 20260923; Wilcoxon *p* = 0.0314; 57 positive and 72 negative pairs, sign test *p* = 0.218). This design had 80% power for d = 0.247.

Within the 281 published pairs, the difference comes mainly from the 56 pairs whose state partner was tested (mean −0.748, nominal *p* = 0.00148), not from the 225 pairs with an untested partner (−0.145, *p* = 0.233). Genes and pairs are the units in these analyses. The sign-flip null does not account for genes that share loci, and none of these *p* values is calibrated. The 1,000-draw matched-pair label control is centered near its analytical expectation of zero; it does not validate the relabel-and-rematch null. These comparisons do not establish a formal biological separation between the two gene classes.

### Protein-response comparison

This comparison asks whether the liver protein response to fibrosis stage differs between the two classes. It contains 73 pairs, with a mean residualized difference of −0.015601, a paired effect d of −0.08156 and a minimum detectable effect of 0.3279 at 80% power. We detected no difference (post hoc paired sign-flip nominal *p* = 0.487465). This is not evidence of equivalence. The relabel-and-rematch *p* and *q* (both 0.459654) are invalid diagnostics under the same null limitation. The unresidualized protein comparison is circular: class selection uses the RNA effect, and protein abundance tracks that effect. It therefore provides no independent biological evidence.

### Expression observability of the two classes

Among genes with a single-cell hepatocyte-specificity value (399 of 400 genetic-only and 1,124 of 1,129 state-only genes), genetic-only genes have mean bulk AveExpr 4.2115, against 0.4164 for state-only genes, a difference of 3.7951 log2 units. Of the 399 genetic-only genes, 281 (70.4%) obtain an expression-matched partner.

## Supplementary Note 4. Fine-mapped coding and splice posterior by trait class and ancestry

Averaged over 1-Mb regions and pooled across ancestries, credible sets for direct MASLD and liver-fat traits carry 23.34% of their posterior probability in protein-altering or canonical-splice variants, against 3.87% for the liver enzymes ALT, AST and GGT. The difference is 0.1947 (95% physical-cluster bootstrap interval 0.0692 to 0.3399).

The pooled contrast mixes trait class with GWAS ancestry and fine-mapping coverage. Within European studies, the fractions are 29.88% and 6.92%, and within admixed American (AMR) studies, 6.73% and 4.54%. No stratum supports a calibrated trait-class call. The European estimate failed the coverage check (the UKBB enzyme GWAS were fine-mapped at 8 of 767 lead windows) and the power check. The European, AMR and combined estimates all failed their random-split calibration controls. The pooled direct-trait excess comes mainly from the three PDFF GWAS and 5 of the 19 direct clusters. The ancestry-stratified estimates used 10,000 shared draws over 81 pooled clusters (seed 20260923). The 57 robustness comparisons gave ratios from 4.37 to 8.28.

In a sensitivity analysis, we used all 6,467 credible sets, including those that failed the PIP-sum criterion (seed 20260817). Matching on standardized lead-variant effect gave a direct-to-enzyme ratio of 5.78 (95% interval 1.99 to 48.05; 27 matched enzyme loci). Direct credible sets come from smaller GWAS (median effective n 21,893 against 160,000, over the 208 direct and 6,259 enzyme credible sets in this analysis) and carry larger standardized effects (median 0.0368 against 0.0183).

To group credible sets into 1-Mb regions, we placed each set at the positions of its member variants in the fine-mapping output. For 715 of the 6,415 credible sets, the position in the locus identifier is too small by a factor of 10 or more and would place the set in the wrong region. All Figure 2D estimates use the member-variant positions.

Separately, the strongest evidence for 423 of the 436 colocalized genes comes from a liver-enzyme trait and for 13 from a direct disease or liver-fat trait. In addition, 22 genes have a direct-trait multi-signal SuSiE colocalization (PP.H4 above 0.5). Multi-signal colocalization tested far more enzyme than direct-trait gene-by-GWAS pairs (12,121 against 771).

## Supplementary Note 5. Signed agreement between ALT TWAS and the disease-versus-control contrast

To test whether genes with higher genetically predicted liver expression for ALT also differ in the same direction between disease and control liver, we compared ALT TWAS z-scores from OTTERS [ref: OTTERS TWAS method] with the disease-versus-control t statistic across 14,703 genes. The TWAS z-scores come from a single OTTERS analysis that was not recomputed with the other genetic analyses in this paper, so this comparison is secondary to the colocalization results.

Without composition adjustment, genome-wide rank agreement between TWAS z and the disease-versus-control t was bounded inside ±0.03 (expression-stratified Spearman rho 0.0107, 95% CI −0.0059 to 0.0272; within-chromosome shift *p* = 0.20). After adjustment for hepatocyte fraction, rho was 0.0167 (95% CI 0.0001 to 0.0332). This interval is not inside ±0.03, and we state no direction. Among the 1,527 genes with the strongest predicted-expression associations (|z| > 3), rho was 0.0421 (95% CI −0.0075 to 0.0916), and 0.514 of genes had the same sign, against 0.499 expected. The unadjusted bound excludes sign agreement of about 0.60 or more among those 1,527 genes. It does not exclude full agreement among the 103 of them that are canonical disease-state genes.

## Supplementary Note 6. Re-annotation of noncoding context by genomic position

Variants absent from the precomputed annotation table received no noncoding context call. Such variants were rare in European studies and common in non-European studies (0.7% of gated credible-set member rows without a call in EUR, 98.2% in AMR). We therefore re-annotated promoter and open-chromatin context by GRCh38 position, with all 16,324 sites lifting over 1:1, and ABC-enhancer context by GRCh37 position. After re-annotation, every noncoding posterior probability had a context call.

## Supplementary Note 7. Sensitivity analyses for the bulk disease-state layer

This note gives the secondary sensitivity statistics for Figure 3. Each section names the main-text result it qualifies. Units are genetic-identity clusters, not verified independent donors.

### GSE126848 per-cohort model

The GSE126848 model contains 23 control and 7 disease clusters, all male, so it omits the constant sex term. Native fragment counts and an independent refit with separate group means reproduce all 6,147 calls. A retained disease library is deeper than its paired control in 5 of 161 control-disease pairs (exact Wilcoxon *p* = 1.9e-5). Each of the 437 calls with zero reads in all 30 libraries has log2FC 0.554. This value equals the precision-weighted control-minus-disease difference in mean log2 normalized library size. Adding log2 library size as a covariate keeps 2,025 of the 6,147 calls (2,678 calls in total) when library size is counted before normalization, and 3,056 (3,649 in total) when it is counted after. Neither refit keeps any of the 437. A cohort-specific expression filter gives 3,440 calls on 15,592 genes.

### Robustness of the 1,854 pooled calls

In the spike-in calibration, 99.8% of supported one-cohort effects (2.0 log2FC) lost their call in at least one leave-one-cohort-out refit. A uniform 1.0 log2FC effect kept its call in all five refits in 90.6% of cases, and losses concentrated in calls near the 0.5 log2FC threshold.

With sample weights fixed at their primary values, adding within-cohort standardized log library size keeps 1,760 of the 1,854 calls when library size is the total count on the tested genes, and 1,648 when it is taken after RLE normalization. Adding top-100-gene share or expression entropy instead keeps 1,738 or 1,759. Adjusting for log library size, top-100 share and detected genes together keeps 1,563 (84.3%) or 1,569 (84.6%) calls under the two library-size definitions, and up calls fall from 1,333 to 1,085 or 1,100.

### Genetic ancestry sensitivity

Adding four ancestry PCs changed 233 calls, against a median of 120 (95th percentile 161) when the same PCs were permuted within cohort and group. Restricting to EUR-called units changed 570 calls (101 lost, 469 gained), against a median of 409.5 (95th percentile 501) for random subsets of the same cohort and group sizes. Both Holm-adjusted *p* values are 0.0199 over the two tests. Across the 1,854 calls, the slope of the with-PC log2FC on the primary log2FC is 0.984 (descriptive). The ancestry-sensitive rule flagged 64 genes; the permutation reference gives a median of 6 expected false flags for this rule (95th percentile 17).

### Transfer of grouped fibrosis-stage effects to GSE213621

The discovery-agreement ratio divides the GSE213621 slope on effects from two discovery cohorts by the held-out discovery cohort's slope on the same effects, with each slope averaged over the three held-out cohorts. It is 0.793 (95% interval 0.488 to 1.098). The F2 versus F0–F1 contrast was untestable, because the minimum detectable slope under the permutation null was 0.34, above the prespecified limit of 0.30. Adding GEO accession quartile as a covariate gives b = 0.719, so the F3–F4 agreement is not attributed to submission order. With estimated cell-type composition added to both fits, b is 0.431, above half of the primary slope. This adjusted result depends on the deconvolution input, and no gene is named from an adjusted fit. In women only (161 units), b is 0.742, against 0.708 in the primary fit.

### Differences at earlier fibrosis steps for GNMT, MAT1A and CYP2C19

GNMT and MAT1A have lower abundance at F1 than at F0 under all three F0 references (GNMT q = 0.00913 / 0.0426 / 0.0430; MAT1A q = 2.45e-4 / 9.5e-4 / 1.6e-3), but not without GSE130970 (q = 0.207 and 0.0892). CYP2C19 has lower abundance at F2 than at F1 (q = 0.0477). Its q is 0.052 with a control-status term, and 0.219 and 0.499 without GSE135251 and GSE162694. GNMT's F1-to-F2 q is 0.0562. We do not use these earlier-step differences to place where the lower abundance begins.

### Protein covariance of GNMT, MAT1A and CYP2C19

Across the 12 pair-by-adjustment tests in the same 46 participants, BH q is 0.00972 to 0.0232 without histology adjustment and 0.0232 to 0.0417 after adjustment for fibrosis stage. After adjustment for NAS, q is 0.136 to 0.214, and after adjustment for the four histology components it is 0.287 to 0.301. We report correlations and pointwise participant-bootstrap intervals for all adjustments, including a prespecified broader metabolic-abundance sensitivity (Figure S3 [panel letter to be assigned]).

## Supplementary Note 8. Stage and lineage sensitivity analyses for the single-cell programs

Units in these analyses are biological donors.

### Equal-variance stage sensitivity

The 62-donor equal-variance fit adds two GSE189600 donors without ambient correction. It gave Ductular injury q = 0.03763 (HC3 q = 0.1558). The call does not hold with the within-dataset permutation null (q = 0.123), without those two donors (q = 0.0985) or with a 50-cell lineage gate (q = 0.0765) (Figure S4 [panel letter to be assigned]).

### Ambient-RNA audit

The Stromal ECM paired change (nominal *p* = 0.0064) has the smallest ordinary *p* among 106 computable paired changes. Its HC3 *p* is 0.024, eighth by HC3 *p*.

### Content-lineage contrasts

Both content-lineage results are reported with cohort-specific intervals, partial correction coverage, and full-raw, 20-cell and annotation-confidence sensitivities (Figure S4 [panel letter to be assigned]). A transfer of the member genes of the two programs to bulk tissue is also reported as a supplementary analysis.

## Supplementary Note 9. Matched-background comparison for Hotspot program associations in bulk liver

In the 364-unit arm of Figure 4E with unweighted scores, 1 of 114 programs reaches fibrosis BH q < 0.05 (228-test family), against a mean of 1.24 for random gene sets matched on size and expression (empirical *p* = 0.473, 200 draws). This unweighted family differs from the weighted 224-test family of Figure 4E. No program reaches q < 0.05 for NAS, against a mean of 0.19 (*p* = 1). In the 457-record sensitivity arm, the counts are 20 against 26.9 (*p* = 0.871) for fibrosis and 2 against 2.57 (*p* = 0.687) for NAS. The programs can be useful for organizing results even though no more of them reach q < 0.05 for fibrosis or NAS than matched random gene sets do.

## Supplementary Note 10. Program-pair covariance in bulk liver and across Visium spots

### Scope

We report program-pair covariance in bulk liver as a discovery-cohort result, with replication in GSE281797 and GSE225740 reported separately for each cohort. Neither replication cohort was checked by genotype for participant overlap with the discovery cohorts. The composition-adjusted analysis in GSE281797 failed its prespecified acceptance rule and was not interpreted. Bulk program scores mix cell abundance with program activity, so these data do not show coordination between cell types. The covariance says nothing about timing or cause.

### GSE213621 test

The same edge sign-concordance statistic could not be read in GSE213621. The prespecified test used 258 GSE213621 disease units (genetic-identity clusters) and 2,328 pairs rediscovered on 371 identity-qualified discovery units, and it was untestable. Its program-group permutation control permutes GSE213621 units within program gene groups, which keeps covariance within each program and removes it between programs. This control still gave a positive mean excess over the non-coordinated baseline of 0.025 (95th percentile 0.081; to four decimals, 0.0251 and 0.0809; prespecified limits 0.02 and 0.05). The negative control therefore failed, and the run stopped before any GSE213621 outcome statistic was computed.

### Replication in GSE281797 and GSE225740

With 3 residual components removed, the excess over random placement was 0.187 (0.118 to 0.188) in GSE281797 and 0.241 (0.186 to 0.242) in GSE225740. These are percentile intervals from 1,000 participant-bootstrap draws within fibrosis stage. Most bootstrap excesses fall below the full-sample value (bootstrap medians 0.162 and 0.217; 4% and 3% of draws exceed the point estimate), so each point estimate lies near the upper end of its interval. Counts were normalized within each cohort; applied to GSE193066, this choice changed the fraction from 0.897 to 0.918. In GSE225740, the result held after adding three estimated composition components, log library size and the top-100-gene count share to the design. The GSE281797 composition estimate failed its prespecified quality rule and was not read. These numbers are not pooled with any GSE213621 or GSE193066 value.

The same pairs were also positive above the random-placement null among participants staged F0-F1 at biopsy. With 1 residual component removed, the excess was 0.166 in 66 GSE281797 participants and 0.304 in 53 GSE225740 participants, whose identities cannot be checked. The combined *p* was 0.000999 for each, and both strata were rejected under Holm correction over the two. This comparison is cross-sectional. In GSE281797, 29 of the 66 F0-F1 participants carry the source diagnosis "No Pathology", so the stratum result cannot be attributed to MASLD remodeling. Without them, the excess was 0.145 (37 participants). Fibrosis stage is already in the model, so the stratum result shows covariance within a narrow stage range. It does not rule out confounding by stage and says nothing about when the covariance arises.

### Discovery-set diagnostics and the GSE193066 comparison

A separate analysis, on a bulk sample set other than the 774 units, reported 2,426 of 6,195 program pairs beyond saturated histology adjustment, modularity 0.358 (against 0.138 for rank-matched and 0.089 for rewired nulls) and a 15-program fibroblast-and-cholangiocyte community. In its comparison in GSE193066, with the model fixed in advance, 0.8974 of coordinated pairs kept a positive sign (89.7% edge-sign concordance), against 0.7095 (70.9%) sign agreement for non-coordinated pairs. Modularity was 0.252, and community density was 0.876 against 0.681. The GSE193066 comparison uses the same excess over non-coordinated pairs as the GSE213621 test. This comparison used 105 source-declared GSE193066 participants, each at the first biopsy only, and their identities are not verified. These values are diagnostics only. They are not independent-cohort replication, and GSE193066 is excluded from primary replication, repeat-biopsy ICC and derived ceiling claims.

### Covariance across Visium spots

Across Visium spots, the 160 gene-disjoint coordinated cross-lineage pairs showed an excess of 0.0257 Fisher z, measured over a label null stratified by lineage pair. The excess remained after adjustment for lobular position, for the first spatial component and for program-specific effects.

## Supplementary Note 11. Equal-cohort pooled heatmap of the 43 molecular systems

The pooled heatmap of the 43 systems balances collections within each cohort and then gives each cohort one-fifth weight (Figure S4U). It shows the profile shared across cohorts without letting the largest cohort dominate. The adjacent panel shows the cohort-specific profiles, including discordant ones. Both panels use the same descriptive pooled late-minus-early row order. This order makes increasing and decreasing shared profiles easy to see, but it is not a support or timing statistic.

## Supplementary Note 12. Retrospective transfer of the continuum loading

GSE174478 and GSE240729 had already been evaluated, so they serve as retrospective transfer checks. These values come from a transfer of the 139-gene continuum loading built on an 844-record bulk set, and they were not recomputed on the 774 units. Sample-level correlations between the continuum and fibrosis were 0.510 (93 samples) and 0.710 (66 samples), with matched-background Holm *p* = 0.117 and 0.064 (1,000 matched random-signature draws), neither below 0.05. A positive association with histology therefore does not show that the loading performs better than matched random signatures. GSE240729 also uses FFPE tissue.

## Supplementary Note 13. RNA-protein agreement in Figure 5B

Across the 3,054 measured genes of Figure 5B, Spearman *rho* between the RNA and protein changes was 0.458, and RNA and protein changed in the same direction for 2,013 genes (65.9%). HLA-DQA1 (Figure 5C) is quantified with the primary-assembly annotation ENSG00000196735.14.

## Supplementary Note 14. Figure 7 severity instrument: fit and scores against task-native baselines

The released linear severity instrument was fitted to 521 sample records from five cohorts. This fit is separate from the 774-unit differential-expression fit of the Resource bulk analysis. All 521 records passed through the RNA-based genetic-identity screen. It found eight within-cohort duplicate pairs (KING kinship > 0.354; three in GSE130970 and five in GSE135251) and no duplicate across cohorts. The Figure 7B fibrosis comparison holds out one cohort at a time, so no duplicated library is in both training and test data. The duplicates remain within the sample bootstrap, which therefore slightly overstates the number of independent people.

The severity instrument scored 0.6631 against 0.6094 for the in-cohort ridge baseline (out-of-cohort Spearman correlation with Kleiner stage in 109 GSE268273 participants). GSE268273 was also used to choose the expression transform. The paired difference, +0.054, has a 95% interval of −0.025 to +0.137, which includes zero. The cis head predicts each H3K27ac region with a ridge regression that adds, to a global expression term, the expression of up to 10 genes whose transcription start sites lie within 50 or 100 kb of the region; window and gene number were chosen inside each training fold. It scored 0.1426 against 0.0480 for a ridge on 81 global expression components without cis genes (out-of-fold skill in 99 GSE267145 participants). This comparator is weaker than the rank-4 reduced-rank regression of Figure 7B (skill 0.122 in the same 99 participants). In a post hoc comparison on the same folds and targets, the cis head exceeded a reduced-rank regression with rank chosen in each fold (skill 0.119) by 0.023 (95% participant bootstrap 0.015 to 0.031). This is an upper bound, because the rank search reached its limit in 3 of 5 folds.

## Supplementary Note 15. Power of the comparison with Kumak et al.

Planted shifts equal to an odds ratio of about 1.74 or 0.48 in the NAS list were detected in 82% and 65% of simulations, and shifts of about 1.75 or 0.48 in the fibrosis list in 77% and 62%.
