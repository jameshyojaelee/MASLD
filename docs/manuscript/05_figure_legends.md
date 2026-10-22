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

**(D)** Fine-mapping posterior in four exclusive consequence classes, averaged over physical clusters (lead positions chained within 1 Mb) for direct MASLD/liver-fat traits (164 credible sets, 16 clusters) and liver-enzyme traits (3,943 sets, 115 clusters); only credible sets whose lead variant has p < 5 × 10⁻⁸ are included. Inset: direct minus enzyme protein-altering plus splice probability, 23.0 percentage points (95% cluster bootstrap 8.1 to 39.6). Ancestry differs between arms, so this is not a trait-class contrast.

**(D, ancestry rows)** The same difference by GWAS ancestry. Thick black bars, 95%; thin gray bars, 98.3% cluster-bootstrap intervals (Bonferroni over EUR+AMR, EUR and AMR). Pooled and all-non-EUR rows are descriptive. No member receives a trait-class call.

**(E)** Ancestry of the GWAS supporting each of the 436 multi-signal genes (SuSiE PP.H4 > 0.5); a per-gene, not per-locus, partition. Of 340 EUR-only genes, 127 were tested without support in non-EUR studies and 213 were untested, which is not the same as unsupported. Black strip, non-EUR ABF PP.H4 > 0.5. The eQTL panel is European, so non-EUR GWAS are not ancestry-matched.

**(F)** All 13 genes with PP.H4 > 0.5 (SuSiE or ABF) in four of five ancestry groups. Circle area, maximum PP.H4 over each group's studies. Dark blue, SuSiE; light blue, ABF only; gray, neither above 0.5. Asterisk, four groups by SuSiE alone. These partitions do not establish biological differences between ancestry groups.

**(G)** RORA in MVP NAFLD EUR (N = 432,709; 29,342 cases). Lead p = 7.6e-8; no fine-mapping record, so no PIP track. SuSiE PP.H4 = 0.983.

**(H)** FABP1 in UKBB ALT (N = 343,850). Lead 2:88423510 (PIP 0.57) and COLOC top variant 2:88424066 (PIP 0.43) are hollow diamonds. SuSiE PP.H4 = 0.998. FABP1 is not differentially expressed in bulk liver (log2FC −0.04, padj 0.51).

**(I)** Key for (G) and (H). Tracks: GWAS p colored by r² to the lead (EUR LD), SuSiE PIP, Broadaway liver eQTL (N = 1,183) [ref: Broadaway et al., liver eQTL] and genes (500 kb, hg19). Gene fill is bulk log2FC, MASLD versus control, over 774 genetic-identity clusters (not verified donors), capped at ±0.3; gray, padj ≥ 0.05; outlined white, outside the bulk gene universe.

See also Figure S2.

## Figure S2

**Figure S2. Fine-mapping and colocalization detail, related to Figure 2.** Panels E and N were withdrawn; the remaining panels keep their letters.

**(A)** The 35 prespecified Tier-1/2 GWAS strata, in 33 rows colored by GWAS ancestry (Pan-UKBB Central/South Asian labeled SAS). Left, loci with a SuSiE or CARMA credible set (log scale); the MVP European ALT, AST and NAFLD strata had no standalone fine-mapping. Right, genes with multi-signal (SuSiE-coloc) or single-signal (coloc.abf) PP.H4 > 0.5. A gene colocalized in several strata counts in each row; the union over the 35 strata is 996 genes.

**(B)** Independent fine-mapped loci. Credible-set windows (1,665 in the 32 fine-mapped strata) are padded by ±250 kb and merged across strata, giving 265 loci (182 in one stratum, 83 in two or more). Method bars are split by best variant PIP. SuSiE and CARMA are counted on the merged loci; SuSiEx and meSuSiE on their own MVP cross-ancestry loci, a different locus set. Without padding there are 289 loci.

**(C)** The 11 merged loci of (B) with a result from all four methods (cross-ancestry loci matched within 500 kb), counted by best PIP. SuSiE and CARMA values are maxima over single-ancestry strata; SuSiEx and meSuSiE values come from one joint MVP fit, so the comparison is not like-for-like. One locus is 9 percentage points, and the panel does not support a difference between methods.

**(D)** EFHD1 locus for AST on chromosome 2 (hg19). Tracks show European (UKBB) and East Asian (BBJ) GWAS association colored by r² to the shared lead, BBJ SuSiE PIP (filled, credible-set member; UKBB AST has no credible set at this locus in the untargeted fine-mapping run, although the colocalization-targeted run has one) and the EFHD1 liver eQTL. Genes are filled by bulk log2FC when they have a canonical call (EFHD1 −1.10) and gray otherwise.

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

**(Q)** 17 genes from a published liver single-cell chromatin study (Elison et al., 2025), grouped by the source's layer. Columns give the maximum eQTL PP.H4 (either method) over the 15 direct MASLD/MASH/liver-fat strata and over the 20 liver-enzyme strata, the published finding and cell type, and the bulk call with BH padj. An asterisk marks PP.H4 > 0.5 from coloc.abf only (multi-signal PP.H4 at most 0.5 or not run).

**(R)** 27 genes whose locus carries an MPRA differential-activity variant in a published liver MPRA (Zhu et al., 2026). Columns as in (Q), plus the strongest allelic MPRA log2FC in HepG2 or LX-2 cells and every map in this study with support. NCAN and DLG5 are outside the 24,115-gene bulk universe and are not evaluable, which differs from no call.

**(S)** Coding plus canonical-splice share of fine-mapped posterior (Figure 2D; credible sets at leads with p < 5 × 10⁻⁸), split within the direct class. PDFF (three European GWAS) and NAFLD plus NASH are compared with all enzyme traits and with each other; PDFF is also compared with European enzyme sets only, because all PDFF sets are European. Bars are 95% physical-cluster bootstrap intervals (10,000 draws; 20,000 with subclasses resampled independently for PDFF minus NAFLD plus NASH).

**(T)** Phenotype provenance of the 35 strata: liver enzyme (ALT, AST or GGT; 20 strata, 5 European), direct MASLD/MASH diagnosis (12; 9 European) and MRI-PDFF or histologic steatosis (3, all European). Non-European GWAS are paired with the European eQTL panel, which limits colocalization.

**(U)** SuSiE-coloc genes (PP.H4 > 0.5) per trait-ancestry cell, and those with a lead within ±250 kb of a fine-mapped window from the same study; the Windows column counts fine-mapped windows. Totals are 436 genes, 100 with such a lead, and 1,665 windows. European enzyme cells have few such genes because MVP European strata were not fine-mapped and UKBB strata only at few loci.

**(V)** Regulatory context of noncoding fine-mapped posterior after re-annotation by position (Supplementary Note 6), as a share of noncoding posterior summed over all credible sets of a class at leads with p < 5 × 10⁻⁸. Contexts overlap and are not summed. Pooled values mix ancestries. Within AMR, no context difference received a call, because every calibration control failed.

**(W)** Target-gene biotype of the 436 primary SuSiE genes: 424 protein-coding and 12 lncRNA. Targets are gene-level and do not resolve transcript isoforms.

Bulk calls use the 774-unit disease-versus-control fit (genetic-identity clusters, not verified donors). A canonical call is BH padj < 0.05 over 24,115 genes and |log2FC| > 0.5.

## Figure 3

**Figure 3. Bulk liver RNA from five cohorts defines disease- and fibrosis stage-associated genes.**

Units are 774 genetic-identity clusters (138 control, 636 disease), not verified donors. The pooled fit tests 24,115 genes and gives 1,854 canonical calls (BH padj < 0.05 and |log2FC| > 0.5).

**(A)** Metadata recorded by each cohort (tiles) and units (bars, disease and control). GSE162694 records NAS for 117 and fibrosis for 109 of 140 units; the 31 without fibrosis are controls.

**(B)** Units by exact fibrosis stage and NAS (364, controls included). The top row counts all 387 units with an exact stage (F0/F1/F2/F3/F4 = 103/102/88/69/25), used in (E) and (F). GSE213621 records only fibrosis bins.

**(C)** Per-cohort calls (BH FDR < 0.05 and |log2FC| > 0.5). Four cohorts use the fixed 24,115-gene family. GSE126848 calls come from a refit on the 15,592 genes that pass its own expression filter, with BH over those 15,592 genes. Top, calls by direction. Middle, genes per overlap class in the 12,180-gene union; darker segments are genes with a GSE126848 call. Bottom, mixed-sign share of multi-cohort genes: 1,211 of 3,839 (31.5%) with GSE126848 and 346 of 2,565 (13.5%) without it. Shared membership is not replication. In GSE126848 (7 disease, 23 control clusters), disease libraries have about half the depth of controls, and GSE135251 calls cannot be matched for genetic ancestry (6 of 10 controls of African genetic similarity).

**(D)** Group-blind PCA of 113 strict controls and 279 advanced-disease clusters (NAS ≥ 5 or fibrosis ≥ F3) in four cohorts, adjusted for cohort and inferred sex. Open gray circles, controls; filled circles, exact fibrosis stage; open triangles, 77 GSE213621 units in the F3–F4 bin. Top, PC1 density by entry rule. These units overlap the pooled fit, so this is a sensitivity analysis, not independent replication.

**(E)** Top, genes differing from F0 at each stage (BH q < 0.05 within each 24,115-gene family, no fold-change filter; 1,338 / 4,266 / 6,190 / 7,043 for F1–F4), with cohort and inferred-sex terms. Lighter inner bars, genes also with |log2FC| > 0.5. F0 holds 14 source controls and 89 disease biopsies. The F1 count rests mainly on GSE130970; without it, 14 genes remain. Bottom, F0-centered activity of five k6 NMF axes (continuous, not patient classes; k6-P6 fails the fixed reconstruction rule and is not shown) with 95% intervals (estimate ± 1.96 SE). F0-reference and cohort-omission views are in Figure S3D–F.

**(F)** Fifty-four genes of the 436-gene primary SuSiE set with BH FDR < 0.05 in at least one adjacent fibrosis contrast (205 / 190 / 157 / 94 clusters). Color, higher-minus-lower stage log2FC, clipped at ±0.688 (triangles). Filled dots, BH q < 0.0125 (5% split across four contrasts; 15 genes); open dots, 0.0125 ≤ q < 0.05. FDR is controlled within each contrast, not over their union. No F0→F1 mark keeps BH q < 0.05 without GSE130970, and no F1→F2 mark without GSE162694. Right, PP.H4 and GWAS trait of each gene's strongest tier-1/2 driver. Selection is conditioned on stage significance; no enrichment or causal claim is made.

**(G)** Disease-state |log2FC| against unsigned SuSiE COLOC PP.H4 for 3,769 genes with a SuSiE test. Dotted lines mark the class rules. The lower strip expands PP.H4 < 0.05, where 2,939 Neither genes are binned. Open circles, seven SuSiE-untestable genes, drawn at their best finite PP.H4 from another tier-1/2 study and not counted as Neither. States precede the shared-posterior check of carried-forward negatives (Supplementary Note 3). Counts use the native RNA gene universe, not the Gene Catalog universe (Figure 6). Fisher and Spearman tests are nominal and unstratified; the prespecified stratified test is in Figure 1B. Inset, lncRNA versus protein-coding: odds of a call over 21,501 genes, unadjusted (open) and adjusted for expression and detection (filled), 95% CI; share of calls up; conditional odds of up in 520 covariate-matched pairs (Wald 95% CI).

See also Figure S3.

## Figure S3

**Figure S3. Bulk transcriptomic sensitivities and continuum analyses, related to Figure 3.**

Except in (H) and (X), units are drawn from the 774 genetic-identity clusters, which are not verified independent donors. GSE213621 records only Control, F0–F1, F2 and F3–F4, so stage displays show it apart from cohorts with exact stages. Category means and overlapping windows are descriptive.

**(A)** Principal component analysis of the 774 clusters (138 control, 636 disease; five cohorts). Input is voom log-CPM of the 2,000 most variable genes after regressing out cohort and inferred sex, without protecting the disease label. PC1 explains 7.9% and PC2 5.2% of variance. Gray open circles are controls and magenta points are disease. The top strip shows PC1 density per group.

**(B)** Log2FC for advanced disease (n = 279) versus strict control (n = 113) in four cohorts, against log2FC from the five-cohort disease-versus-control fit (774 units). Points are genes called in either contrast; dashed line, y = x. Over all 24,115 genes, Spearman rho = 0.953. Of 1,854 canonical calls, 1,755 are also called in the endpoint contrast (recovered, magenta) and 99 are not (subthreshold, pink). A further 1,195 are called only in the endpoint contrast (gray). No gene is called in opposite directions.

**(C)** Direction agreement of the 1,854 pooled calls in each cohort ranges from 78.3% to 100% (GSE126848, 84.36%). Leave-one-cohort-out fits on the fixed 24,115 genes give effect Spearman correlations of 0.870–0.995 and call retention of 70.7%–96.6%. Full refits that repeat filtering after cohort removal keep 66.7%–96.8% of calls as same-sign calls; among calls still tested, 99.7%–100% keep their sign.

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

**(C)** Ordinal stage β for the 117 programs, ranked within lineage, from generalized least squares with per-dataset residual variance in 60 ambient-correctable donors (GSE202379, GSE244832; donors per lineage above each block). 17 of 53 T-cell donors have fewer than 20 T cells; without them (36 donors), no T-cell program is selected (smallest q over 117 = 0.587). Gray bars, β ± standard error times the within-dataset permutation critical |t| at 0.05/117; black bars, nominal 95% confidence intervals. No program passes the 117-program permutation BH and stability rule (smallest q = 0.323).

**(D)** Discovery against content-associated lineage for the 117 programs (donor-level same-atlas scores). Eighty match (shaded diagonal) and 37 do not. Labels mark the 7 calls that keep the same lineage in fewer than 95% of 1,000 donor bootstrap draws. Content association does not establish transcript origin or cell autonomy.

**(E)** Fibrosis effect (cohort-standardized program-score SD per recorded stage-code unit) for 112 eligible programs, ranked, in 364 genetic-identity clusters; NAS effects in Figure S4W. Filled points, random-effects meta-estimates colored by discovery lineage; gray bars, Knapp-Hartung 95% confidence intervals (t, 2 df); open marks, per-cohort estimates (circle GSE130970, n = 71; triangle GSE135251, n = 207; square GSE162694, n = 86). No association meets BH q < 0.05 over 224 tests (112 programs x fibrosis and NAS; smallest q 0.079 and 0.159), with or without composition adjustment. The intervals are not adjusted for the 224 tests, and with three cohorts their width depends mainly on between-cohort agreement (Q), so 34 of the 112 fibrosis intervals (and 12 NAS intervals, Figure S4W) exclude 0 without a BH call. With z-based REML inference instead, 37 fibrosis and 14 NAS associations reach q < 0.05.

**(F)** Cohort-standardized Stromal ECM and Ductular injury scores in nine overlapping fixed-projection windows (774 genetic-identity clusters, five cohorts). Lines join window means; bands are ± 1 standard error; clusters per window in key. Dashed gray lines, cohorts overlapping continuum construction; solid lines, the two independent evaluation cohorts. The Stromal ECM score excludes 9 continuum-signature genes carrying 19.5% of its weight. Descriptive; window means do not imply order in time. Adjusted estimates in Figure S4V.

**(G)** Stage-adjusted continuum coefficient for each of 5,435 testable nodes (gray; pathways and Hotspot programs) in 43 fixed, signature-excluded systems, estimated in 367 stage-complete disease units (109 GSE162694, 258 GSE213621). Rows, systems by size (node count in brackets). Open circles, system medians with 95% bootstrap intervals (disease units resampled within cohort, B = 1,000). Four nodes are untestable in the 24,115-gene data, because fewer than 80% of their original members remain observed after signature exclusion; they are counted at right.

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

**(F)** Left, illustrative H&E overlay of the adjusted Stromal ECM score (section JBO019, donor H37), chosen outcome-blind; no pathologist annotation. Right, Moran excess over matched genes after adjusting spot scores for all 16 cell2location abundances, in Guilliams et al. and Vu et al. (2026). Small points, donors (Guilliams et al., 4) or arrays (Vu et al., 10); band, matched-null 95% range. Labels give BH q from 9,999 matched-gene draws (two programs per source) and donor or array counts. The gene null holds donors and arrays fixed and draws control genes independently, matched on correlation with hepatocyte abundance, so it is not a donor-level test and is not matched for program co-expression. In a descriptive analysis that was not prespecified, 68 of 114 (Guilliams et al.) and 74 of 114 (Vu et al.) scorable fixed programs exceed a null built the same way, with each program's controls matched on its own discovery lineage (BH q < 0.05 over 114 programs per source); Stromal ECM ranks 42nd and 29th by excess. With controls matched on correlation with fibroblast abundance, Stromal ECM still exceeds the null (q = 0.011 and 0.0008); with controls matched on correlation with cholangiocyte abundance, Ductular injury does not (q = 0.10 and 0.50). Donor and array counts are of units with positive excess over the matched null: Stromal ECM excess is positive in 3 of 4 Guilliams et al. donors (exact one-sided sign-flip p = 0.19; the smallest possible p with 4 donors is 0.0625) and Ductular injury excess in 3 of 4 (p = 0.125). Spot scores use the measured member genes, which carry 55% (Guilliams et al.) and 53% (Vu et al.) of the Stromal ECM L1 weight and 41% of the Ductular injury weight. IGFBP7 carries 2.47% of the full Stromal ECM weight and BICC1 5.44% of the full Ductular injury weight, so the result concerns program scores, not either gene. The two sources used different cell2location reference builds (Methods). Vu et al. has no donor key (open point, source-dependent). The Guilliams et al. donors are surgical samples taken for other conditions, not a MASLD cohort; H37 and H38 are also in the single-cell atlas.

## Figure 6

**Figure 6. The MASLD Gene Catalog keeps eQTL-colocalization and bulk RNA calls separate for 28,233 genes.**

**(A)** Catalog structure. Each gene identity carries a separate eQTL-colocalization call and bulk RNA call, each qualified by testability and provenance, plus one molecular object (protein-coding gene, lncRNA transcript, ambiguous lncRNA locus or untestable nomination; no record is a DNA element) and one proposed discriminating experiment. Protein, chromatin and spatial context keep their assay or program unit and enter no class or route.

**(B)** Genetic call (rows) by bulk RNA call (columns) for 28,233 gene identities. The genetic call is supported when SuSiE-COLOC PP.H4 > 0.5 in any of 35 Tier-1/2 GWAS (15 direct MASLD, 20 liver enzyme). The RNA call is supported when BH q < 0.05 (over 24,115 genes) and |log2FC| > 0.5 in one limma-voom fit of 774 genetic-identity clusters (138 control, 636 disease), which are not verified donors. Blue cells are tested in both domains; gray cells are untested in at least one. No negative-call rule was defined, so there is no tested-negative row or column, and untested does not mean negative. Genetics untestable means no Tier-1/2 COLOC record (9,312) or SuSiE-untestable in one study (7; finite PP.H4 elsewhere, maximum 0.164). The untestable-by-untestable cell is 0 by construction. "Coding PIP ≥ 0.5" counts 31 genes with a fine-mapped coding variant at SuSiE PIP ≥ 0.5. It does not enter the genetic call, so PNPLA3, TM6SF2 and MTARC1 stay unresolved.

**(C)** Four worked entries. Columns give SuSiE and ABF PP.H4 with the top study, any fine-mapped coding variant, bulk RNA log2FC with 95% CI (log2FC ± 1.96 SE; 774 genetic-identity clusters) with BH q, and the Catalog class. The line under each gene gives its gene-level and molecular-object routes. Dashed lines mark |log2FC| = 0.5. ABF PP.H4 is shown for reference and makes no call. *GLP1R* is untested because its RNA was not detected; its genetic call is indeterminate (no Tier-1/2 SuSiE result). *MTARC1* has no Tier-1/2 SuSiE result, so its genetic call is indeterminate. *MTARC1* and *THRB* have indeterminate RNA calls because the change is below the 0.5 floor. The *MTARC1* missense variant (PIP 0.76) does not enter its genetic call. "Both supported" (*HKDC1*) does not compare direction. Genetic support does not establish a causal mechanism.

**(D)** Genes per Catalog class, with the gene-level route printed beside each class: untested 13,551, unresolved 13,127, established-state associated 1,129, genetically anchored 400 and supported on both axes 26 (Catalog class "concordant"). The untested bar splits into no RNA call (4,232, not detected) and no genetic call (9,319). Gray is untested and blue tested. Hatching marks the 708 untested genes with one supported call (698 RNA; 10 genetics, narrower than a line). The 20 unresolved genes with a coding variant at PIP ≥ 0.5 are routed by their eQTL-colocalization class. Of the unresolved, established-state and untested genes, 290, 18 and 25 have ABF PP.H4 > 0.5 but no Tier-1/2 SuSiE result. Routes are proposals, and counts are not a priority ranking.

**(E)** Ways to use the Catalog: view the worked entries, inspect each molecular object with its uncertainty, read source references, and download the gene calls with the molecular-object records. All views draw on the same inputs and give no combined score or rank. The Catalog files are available at [accession to be supplied] (see Data availability).

## Figure 7

**Figure 7. Predictive models are scored against task-specific baselines beside a record of their training exposure.**

**(A)** Evaluation framework. Each of 8 registered tasks has a baseline set in advance. A score whose exposure record is unresolved is uninterpretable, not evidence of skill. The bar counts audit entries by checkpoint finding (104 audits of 106 model entries), before eight later audit revisions. These checkpoint findings are separate from the study-specific exposure in (C) and (F).

**(B)** Gain over each task's own baseline, named with its score on each x axis. Reporter assay (GSE281364, HepG2 contexts only): Bonferroni-corrected intervals across the seven rows (paired block bootstrap); sequence-based reporter baselines were not run. eQTL |z| and sign: held-out LD fold, paired 95% LD-block bootstrap. Fibrosis: samples are the unit because the training pool has no donor key; 95% sample bootstrap. H3K27ac: 95% participant bootstrap. ATAC: donor jackknife.

**(C)** Donor-class-balanced macro F1 for broad cell classes on all seven held-out studies (hollow; 102 donors, 50,000 cells) and on each encoder's clean studies (filled). Gray: HVG+PCA with the same head. Intervals are 95% study-stratified donor-cluster bootstraps. Geneformer V2 has no clean study, so its scores are not evidence of skill. Hollow and filled scores use different study sets and are not compared with each other; (F) compares each encoder with HVG+PCA on the same studies.

**(D)** Signed eQTL: best head's sign AUC per |z| stratum (95% LD-block bootstrap, uncorrected across five strata) beside within-block permutation chance and the sequence-attributable bound, which assumes at most one causal variant per gene and LD block. Severity transfer: released instrument, no refit; Spearman with stage (stage-stratified 95% bootstrap) beside each cohort's stage tie ceiling. GSE213621 units are genetic identities, not verified donors. GSE276114 includes participants with HCC, sampled from peritumoral liver without a per-sample label. LD-block join: Pearson r of each row's effect with its gene's mean effect, own LD block left in or held out. Routing: fitted rule and oracle router against the pre-declared floor.

**(E)** Each statistic is oriented so that positive is the prespecified direction and placed on one scale, (estimate − null)/(threshold − null): 0 is the null, 1 the threshold its verdict uses (minimum detectable effect for severity arms and the zone task; h2 z = 7 for GWAS rows; twice the test's resolution for signed eQTL; the pre-declared floor for routing). Severity-arm intervals (95% stratified bootstrap) show precision only; verdicts use the label permutation null, and an underpowered arm is untestable, not negative. Signed eQTL direction, recorded as tested negative, is drawn untestable because its sequence-attributable bound (open diamond) is below the threshold. Treatment response uses the GSE83452 diet arm (6 resolved vs 16 persisted; permutation p 0.23), which its source analysis records as negative and underpowered.

**(F)** Top: exposure of each held-out study (donors in parentheses) for each encoder. Middle: macro F1 gain over HVG+PCA with the same head, on all seven studies (hollow) and clean studies (filled). Bottom: encoder minus the scVI liver encoder trained here, on the encoder's clean studies. Intervals are 95% donor-cluster bootstraps, uncorrected across rows. In the unweighted study mean, a study with 2 donors and 149 cells gets one seventh of the weight; two of twelve model-by-head blocks change sign when weighted by donors, so no sign claim rests on the unweighted mean.
