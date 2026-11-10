# Figure legends

Genetic counts and dependent overlap, accessibility and Catalog summaries
describe a provisional historical map. Saved eQTL fits for 90 of 436 primary
genes contain duplicate variant identifiers; qualified full refits remain
required. The other 346 genes are not a validated replacement set. RORA and
FABP1 also fail the regularized reference-LD matrix check, so their historical
locus displays are not valid-LD controls.

## Figure 1

**Figure 1. A multimodal human MASLD Resource keeps inherited liver-trait evidence separate from disease-state RNA.**

**(A)** Resource arms, each shown with its biological unit, phenotype provenance, source reuse and assay applicability. Genetics: 35 Tier-1/2 GWAS strata, direct disease traits kept apart from liver-enzyme proxies. Bulk RNA: five cohorts, 774 genetic-identity clusters (138 control, 636 disease), not verified donors. Single cell: 7 datasets, 104 source donors, 102 analyzed; the inferential n differs by contrast. Protein, chromatin and spatial data add context in native units. All arms feed the MASLD Gene Catalog (28,233 gene identities). Reuse of a public source is not independent validation.

**(B)** Overlap of liver-eQTL-colocalized genes (SuSiE PP.H4 > 0.5) and canonical disease-versus-control DEGs (BH q < 0.05 over 24,115 genes, |log2FC| > 0.5; 774 genetic-identity clusters), among 3,589 genes evaluated after shared-posterior eligibility qualification. Counts: 426 colocalized, 209 DEGs, 26 shared. Gray bars show the exact null within 10 AveExpr bins; darker bars mark its central 95% (14 to 30). Dashed line, expected overlap (21.4); filled point, observed (26); open circle, unstratified expectation (24.8). Lower row: Mantel-Haenszel OR 1.28 (95% CI 0.83 to 1.99; exact stratified two-sided p = 0.28; original test definition, updated eligibility), log axis. Historical expression-matching diagnostics were not rerun in this universe; no direction is inferred. Outside the comparison are 935 canonical DEGs without a multi-signal result, 11 priority-untestable DEGs and 699 without a primary COLOC record (503 lncRNA). They are not genetic negatives.

**(C)** Trait-class support of the 436 colocalized genes (SuSiE PP.H4 > 0.5 in any of 35 primary GWAS): ALT/AST/GGT (20 studies, 424 genes) and NAFLD/NASH/PDFF (15 studies, 22 genes); ten overlap. Other-class segments show support (10 each), finite recorded PP.H4 ≤ 0.5 (46 enzyme-only, 11 direct-only) or no finite recorded SuSiE output (outlined; 368 enzyme-only, 1 direct-only). Recorded availability is unequal: 3,888 genes in enzyme studies, 392 in direct studies and 362 in both, with 56 enzyme-supported and 21 direct-supported among the latter. Availability uses stored r2 output before qualification of carried negatives; it is not current gene-level applicability or a pair-specific eligibility check. No recorded output does not mean negative. Counts are genes.

**(D)** Analysis flow across Figures 1 to 7. The same source-qualification rules (phenotype, biological unit, coverage and applicability, source dependence) apply to every analysis, and each reports a machine-readable evidence state, the unresolved alternative and the next discriminating experiment. Inherited liver-trait evidence (Figure 2), multicohort disease-state remodeling (Figure 3), donor-level programs (Figure 4) and molecular and tissue context (Figure 5) lead to the Gene Catalog (Figure 6) and model evaluation (Figure 7). Arrows show analysis order, not causal links.

## Figure 2

**Figure 2. Liver-trait GWAS colocalize with liver eQTL at 436 genes, mostly through liver-enzyme traits.**

**(A)** Three unlinked counts over 35 Tier-1/2 GWAS strata: studies by ancestry; fine-mapped loci (360 trait-specific, 265 unique) by ancestry of top PIP, with 54 tied loci as a separate segment; and 996 colocalized genes (436 multi-signal, 560 single-signal only), each counted once in the band of its largest PP.H4.

**(B)** Largest multi-signal PP.H4 per gene against SuSiE PIP of the same lead variant in the same GWAS (n = 866 genes; 327 colored, 539 gray outside the display filter). Shape gives lead consequence. Genes whose largest posterior comes from the unfine-mapped MVP European strata have no PIP and are not drawn; absence is not a negative call.

**(C)** Consequence class of one lead COLOC variant per gene (996 genes, 933 variants), first matching hg19 knownGene class in the order coding or splice, UTR, promoter, intron, intergenic; coding includes synonymous. A gene-level summary, not the posterior composition of credible sets.

**(D)** Fine-mapping posterior in four exclusive consequence classes, averaged over physical clusters (lead positions chained within 1 Mb) for direct MASLD/liver-fat traits (164 credible sets, 16 clusters) and liver-enzyme traits (3,943 sets, 115 clusters); only credible sets whose lead variant has p < 5 × 10⁻⁸ are included. Inset: direct minus enzyme protein-altering plus splice probability, 23.0 percentage points (95% cluster bootstrap 8.1 to 39.6). Ancestry differs between arms, so this is not a trait-class contrast.

**(D, ancestry rows)** The same difference by GWAS ancestry. Thick black bars, 95%; thin gray bars, 98.3% cluster-bootstrap intervals (Bonferroni over EUR+AMR, EUR and AMR). Pooled and all-non-EUR rows are descriptive. No member receives a trait-class call.

**(E)** Ancestry of the GWAS supporting each of the 436 multi-signal genes (SuSiE PP.H4 > 0.5), counted per gene, not locus. Of 340 EUR-only genes, 127 have a finite recorded non-EUR SuSiE output at PP.H4 ≤ 0.5 and 213 have none. Availability describes stored r2 output before qualification of carried negatives; it does not establish current pair-specific eligibility, and absent output is not a negative result. Black strip, non-EUR ABF PP.H4 > 0.5. The eQTL panel is European, so non-EUR GWAS are not ancestry-matched.

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

**Figure 3. Bulk liver RNA defines disease- and fibrosis stage-associated genes, with held-out stage replication in SteatoSITE.**

Genetic counts in Figures 1–3 and 6 describe the historical selected map.
Duplicate variant identifiers were found in eQTL fits for 90 of the 436
primary genetic genes. Their genetic calls, overlap counts and downstream
accessibility links remain provisional pending corrected full fits. This
qualification does not alter the tissue RNA effect estimates.

Units are 774 genetic-identity clusters (138 control, 636 disease), not verified donors. The pooled fit tests 24,115 genes and gives 1,854 canonical calls (BH padj < 0.05 and |log2FC| > 0.5).

**(A)** Metadata recorded by each cohort (tiles) and units by group (bars): 774 genetic-identity clusters in the five discovery cohorts and 493 staged SteatoSITE donors (held out, below the gap), 1,267 units in total. Early / mild, F0–F1, or NAFL in GSE126848, which has no stage; advanced / severe, F2–F4, or NASH. Partial / inferred: sex inferred from expression (GSE213621, GSE135251, SteatoSITE), NAS for 117 of 140 GSE162694 units, and fibrosis recorded only as bins (GSE213621 F0–F1, F2, F3–F4; SteatoSITE F0–F1, F2, F3, F4). SteatoSITE records no NAS, age or controls; its 39 source-reported normal livers cannot be identified and are counted in F0–F1.

**(B)** Units by exact fibrosis stage and NAS score (364, controls included); Total row and column are the margins. GSE213621 records only fibrosis bins and SteatoSITE records no NAS, so neither enters the grid. The stage fits in (E) and (F) use the 387 units with an exact stage (F0/F1/F2/F3/F4 = 103/102/88/69/25).

**(C)** Per-cohort calls (BH FDR < 0.05 and |log2FC| > 0.5; 17,347 calls) flowing from each cohort (left) to the number of cohorts calling the gene (right). Four cohorts use the fixed 24,115-gene family. GSE126848 calls come from a refit on the 15,592 genes that pass its own expression filter, with BH over those 15,592 genes. Of the 12,180 called genes, 8,341 are called in one cohort; labels give the genes in each multi-cohort class and the share called in the same direction in every cohort (2 cohorts 65%, 3 cohorts 76%, 4 cohorts 89%, 5 cohorts 85%). Overall, 1,211 of 3,839 multi-cohort genes (31.5%) have mixed signs, and 346 of 2,565 (13.5%) without GSE126848. GSE135251 is one of the cohorts in 325 of those 346 mixed-sign genes. Shared membership is not replication. Within-cohort comparison in Figure S3 [panel letter to be assigned]. In GSE126848 (7 disease, 23 control clusters), disease libraries have about half the depth of controls, and GSE135251 calls cannot be matched for genetic ancestry (6 of 10 controls of African genetic similarity).

**(D)** Group-blind PCA of 113 strict controls and 279 advanced-disease clusters (NAS ≥ 5 or fibrosis ≥ F3) in four cohorts, adjusted for cohort and inferred sex. Open gray circles, controls; filled circles, advanced disease colored by fibrosis stage (0–4). The 77 GSE213621 units recorded only in the F3–F4 bin are colored at 3.5. Top, PC1 density of controls (gray) and advanced disease (magenta). These units overlap the pooled fit, so this is a sensitivity analysis, not independent replication.

**(E)** Left, top, genes differing from F0 at each stage (BH q < 0.05 within each 24,115-gene family, no fold-change filter; 1,338 / 4,266 / 6,190 / 7,043 for F1–F4), with cohort and inferred-sex terms. Lighter inner bars, genes also with |log2FC| > 0.5. F0 holds 14 source controls and 89 disease biopsies. The F1 count rests mainly on GSE130970; without it, 14 genes remain. Left, bottom, F0-centered activity of five k6 NMF axes (continuous, not patient classes; k6-P6 fails the fixed reconstruction rule and is not shown) with 95% intervals (estimate ± 1.96 SE). F0-reference and cohort-omission views are in Figure S3D–F. Right, held-out replication in SteatoSITE (ENA PRJEB58625; archival FFPE total RNA; 493 donors, F0–F1/F2/F3/F4 = 255/74/75/89; F0–F1 is inferred from the source stage code). Points, ordinary least-squares slope of SteatoSITE log2FC (`~ batch + inferred_sex + stage`) on discovery log2FC over discovery BH calls kept in SteatoSITE, with 95% bootstrap intervals (200 draws); text, genes and share with the discovery sign. Filled, testable (minimum detectable slope ≤ 0.30); open, untestable and descriptive. The top three rows are adjacent steps; the bottom two are grouped-stage contrasts. F3–F4 versus F0–F1 replicates (b = 1.462, 95% CI 1.323–1.601; 93.5% same sign, permutation null 52.6%; Holm p = 0.002 over the two grouped contrasts). Discovery for the grouped contrasts is the 371-unit disease-only three-cohort fit; for adjacent steps, the exact-stage fits in (E). Dashed line, slope 1.

**(F)** Fifty-four genes of the 436-gene primary SuSiE set with BH FDR < 0.05 in at least one adjacent fibrosis contrast (205 / 190 / 157 / 94 clusters). Color, higher-minus-lower stage log2FC, clipped at ±0.688 (triangles). Filled dots, BH q < 0.0125 (5% split across four contrasts; 15 genes); open dots, 0.0125 ≤ q < 0.05. FDR is controlled within each contrast, not over their union. No F0→F1 mark keeps BH q < 0.05 without GSE130970, and no F1→F2 mark without GSE162694. Right, PP.H4 and GWAS trait of each gene's strongest tier-1/2 driver. Selection is conditioned on stage significance; no enrichment or causal claim is made.

**(G)** Disease-state log2FC against unsigned SuSiE COLOC PP.H4 for the 3,589 genes evaluated after shared-posterior eligibility qualification under the priority rule. Groups: Neither, 2,980; Transcriptomics only, 183; Genetics only, 400; Both, 26. Canonical RNA calls require BH padj < 0.05 and |log2FC| > 0.5; genetic support requires PP.H4 > 0.5. Without expression stratification, Both contains 26 observed against 24.8 expected (Fisher OR 1.058, nominal p = 0.742). Spearman rho with |RNA log2FC| is 0.00018 (nominal p = 0.991). Figure 1B gives the stratified overlap test. PP.H4 is unsigned, so RNA sign is not a risk-allele direction. Priority-untestable genes are excluded, including genes with a finite result in another study. Counts use the native RNA universe, distinct from the Gene Catalog (Figure 6).

See also Figure S3.

## Figure S3

**Figure S3. Bulk transcriptomic sensitivities and continuum analyses, related to Figure 3.**

Except in (H) and (X), units are drawn from the 774 genetic-identity clusters, which are not verified independent donors. GSE213621 records only Control, F0–F1, F2 and F3–F4, so stage displays show it apart from cohorts with exact stages. Category means and overlapping windows are descriptive.

**(A)** Principal component analysis of the 774 clusters (138 control, 636 disease; five cohorts). Input is voom log-CPM of the 2,000 most variable genes after regressing out cohort and inferred sex, without protecting the disease label. PC1 explains 7.9% and PC2 5.2% of variance. Gray open circles are controls and magenta points are disease. The top strip shows PC1 density per group.

**(B)** Log2FC for advanced disease (n = 279) versus strict control (n = 113) in four cohorts, against log2FC from the five-cohort disease-versus-control fit (774 units). Points are genes called in either contrast; dashed line, y = x. Over all 24,115 genes, Spearman rho = 0.953. Of 1,854 canonical calls, 1,755 are also called in the endpoint contrast (recovered, magenta) and 99 are not (subthreshold, pink). A further 1,195 are called only in the endpoint contrast (gray). No gene is called in opposite directions.

**(C)** Historical fixed-family direction agreement of the 1,854 pooled calls ranges from 78.3% to 100% (GSE126848, 84.36%). This panel predates native GSE126848 filtering: the current comparison has 937/1,056 testable genes with the pooled sign and 798 untestable genes. These within-discovery comparisons are not independent replication. The full-family REML/modified Hartung–Knapp sensitivity gives zero BH calls (minimum q=0.171430; Supplementary Note 19). Leave-one-cohort-out fits on the fixed 24,115 genes give effect Spearman correlations of 0.870–0.995 and call retention of 70.7%–96.6%. Full refits that repeat filtering after cohort removal keep 66.7%–96.8% of calls as same-sign calls; among calls still tested, 99.7%–100% keep their sign.

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

**([panel letter to be assigned])** Called-gene lists across the five Figure 3C cohorts compared with random groups of one cohort. Within one cohort, the 493 SteatoSITE donors were split 50 times into five random groups of about 99, and F3–F4 versus F0–F1 was fitted in each group (`~ batch + inferred_sex + stage`, BH q < 0.05 within group; points, splits; bars, medians). Left, share of genes in the union of calls that are called in one group only: 26.8% within SteatoSITE (2.5–97.5% range 24.8–29.9%) against 68.5% (8,341 of 12,180) across cohorts. Right, share of genes called in two or more groups that are called up in one and down in another: 0.2% against 31.5% (1,211 of 3,839). The SteatoSITE contrast is fibrosis stage, not disease versus control, and its groups share protocol and sequencing batches, so this is a descriptive comparison. It does not identify how much of the Figure 3C disagreement is caused by sampling, protocol or case mix.

**([panel letter to be assigned])** Exonic and intronic change with fibrosis stage in SteatoSITE (493 donors, F3–F4 versus F0–F1, `~ batch + inferred_sex + stage`). Each point is one of 12,010 genes with intronic signal; x, change in exonic reads (log2); y, change in intronic reads, which mark pre-mRNA (log2); dashed line, equal change. Spearman 0.95. Magenta, 114 genes whose exon/intron ratio changes by at least 0.5 log2 (approximately 1.41-fold) with BH q < 0.05 over 12,010 (limma on each donor's log2 exon/intron ratio, adjusted for batch and inferred sex). Blue, named genes; their intronic change is 90–113% of their exonic change. Aligned intronic and exonic effects are consistent with transcription-associated tissue remodeling. Cell mix also moves both read types together, and RNA processing and stability remain unresolved. CYP2C19 and EPCAM have no intronic region free of other genes and are not shown.

**([panel letter to be assigned])** Figure 3E NMF axes by fibrosis stage in the discovery cohorts (blue; 387 units, F0 and F1 combined as reference, `~ dataset + inferred_sex + stage`) and in SteatoSITE (magenta; 493 donors, `~ batch + inferred_sex + stage`). Activity is in standard deviations of each axis score relative to F0–F1; standardization uses all 774 discovery units and the 493 staged SteatoSITE donors, so magnitudes are on different scales. F3–F4 versus F0–F1, one-sided in the discovery direction, 1,000 permutations, Holm over five axes: P1 to P4 p = 0.005 each. P5 is not interpretable: its discovery F3–F4 effect is near zero (−0.01) and its genes depend on total-RNA versus poly(A) library preparation. P6 failed reconstruction in discovery and is not shown.

**([panel letter to be assigned])** Native-cohort transportability of the fixed
1,854 pooled tissue calls. Filled blue points show same-sign effects; open
magenta points additionally meet the cohort's call threshold. Denominators are
the testable fixed calls: 1,854 in four cohorts and 1,056 in GSE126848, with
798 untestable there. These cohorts contributed to discovery, so agreement is
descriptive and does not constitute independent replication.

**([panel letter to be assigned])** Disease-minus-control log2 fold changes
for CYP2C19, GNMT and MAT1A in each cohort and the REML sensitivity. Cohort
intervals are estimate ±1.96 moderated SE; summary intervals use modified
Hartung–Knapp uncertainty with variance factor at least one. These pointwise
intervals are not multiplicity-adjusted: no gene passes BH over the complete
24,115-gene sensitivity (minimum q=0.171430). The examples were specified before
this sensitivity, after their selection in the manuscript.

## Figure 4

**Figure 4. Fixed single-cell programs, lineage content and bulk expression gradients.**

The atlas census, UMAP and atlas-derived estimates shown here use historical
run-level inputs that include repeated physical cells and 41 antibody-tag
library rows. They are not corrected unique-cell results. All 117 program
definitions and weights remain fixed; verification of unique-cell donor effects
is incomplete, and fixed-weight expression projections are not native Hotspot
scores. Independent bulk and external-cohort measurements remain separate.

Program memberships and weights remain fixed. Single-cell inference uses source donor identifiers; bulk inference uses selected genetic-identity clusters, not verified independent donors. All associations are cross-sectional.

**(A)** Seven datasets contribute 275 libraries, 104 source donors and 1,232,318 cells; 273 libraries, 102 donors and 1,232,285 cells enter the analyzed atlas. The 117 programs were discovered in five lineage subsets. October 26 source corrections exclude GSE202379 P98 from disease inference because its control diagnosis conflicts with SAF S1A1F1, retain P30 as a qualified nonsteatotic reference (S0A1F0), and correct GSE244832 D11/D14 to their deposited RNA diagnoses (MASL/MASH). RNA diagnoses do not verify RNA–ATAC pair identities. Histology, diagnosis and etiology remain separate.

**(B)** The unchanged uniform 50,000-cell UMAP sample (seed 42), colored by cell type. Cells are descriptive observations, not independent replicates. Disease-density overlays are omitted because the diagnosis groups contain different dataset mixtures. Colors shared with E denote discovery lineage.

**(C)** All 117 ordinal-stage coefficients, ranked within the five discovery lineages (Hep, hepatocytes; Fib, fibroblasts; Chol, cholangiocytes; Mac, macrophages; T, T cells). Score ~ clinical category + dataset is fitted by REML generalized least squares with per-dataset residual variance; Healthy/Steatosis/Steatohepatitis are coded 0/1/2. At most 59 source-qualified donors enter the primary universe; the two multi-stage datasets supply the slopes (54 hepatocyte, 42 fibroblast, 37 cholangiocyte, 24 macrophage and 52 T-cell donors). Black intervals are nominal 95% CIs; gray intervals use the within-dataset permutation critical |t| at 0.05/117. With 19,999 permutations (seed 20260923), none passes BH over 117 programs plus the unchanged stability rule (minimum q=0.18135). The variance model was chosen after focal slopes were known. Ambient-corrected GLS selects none of 106 computable programs within the same 117-program family (minimum q=0.11115); the 11 corrected T-cell scores are untestable. Nonsignificance does not establish absence.

**(D)** Discovery versus content-associated lineage for all 117 programs: 80 match and 37 differ. Gray cells mark the diagonal. Each small open circle denotes one of seven assignments retained in fewer than 95% of 1,000 joint donor-bootstrap draws; 76 matches and 34 mismatches are stable. Content association does not establish transcript origin. In an exploratory November 2 check, Hotspot without GSE244832 recovered neither focal set (best top-50 overlaps: Stromal ECM 2 genes; Ductular injury 0). Their genes co-vary across 45 GSE202379 donors (Holm p=0.002 for each), but so do 27 of the other 28 hepatocyte programs. Their discovery depends on GSE244832; this covariance does not establish hepatocyte origin.

**(E)** Fibrosis coefficients for 112 eligible weighted programs, ranked, in 364 clusters from GSE130970/GSE135251/GSE162694 (71/207/86). Colored points are random-effects meta-estimates; gray intervals are Knapp-Hartung nominal 95% CIs (t, 2 df); open circle/triangle/square marks show the three cohort estimates, respectively. Units are cohort-standardized score SD per recorded stage-code unit, adjusted for NAS and sex. No fibrosis or NAS association passes BH over 224 tests (minimum q=0.079/0.159); NAS is in Figure S4W. Thirty-four fibrosis intervals exclude zero without a family-level call. Three-cohort Knapp-Hartung inference is low-powered; z-based inference or adding GSE174478 changes the calls. Colors indicate discovery lineage, not content or transcript origin.

**(F)** Ductular injury and Stromal ECM means ±1 SE across nine overlapping 20-percentile windows centered at 10–90%, using the 774-cluster bulk substrate. Dashed gray lines are the three continuum-construction-overlap cohorts; solid purple/orange lines are GSE162694/GSE213621. Profiles are descriptive, not temporal. In the 367 disease units, adjusted coefficients are 0.707 for Ductular injury and 0.760 for Stromal ECM, in score SD per continuum SD. Exploratory November 2 expression-matched random axes gave medians 0.538/0.612 and 95th percentiles 0.736/0.804, respectively; neither focal association exceeded this reference after Holm over four tests (p=0.43–0.46). Random axes called a median 93/117 programs against 94 observed. Thus the associations do not establish continuum-signature specificity. This control tested programs, not the 43 systems.

**(G)** All 43 fixed systems, continued across two columns in registry size order (S01–S43; full names in source data). Gray points show 5,435 testable node coefficients; open circles and bars show system medians and nominal 95% donor-bootstrap intervals (1,000 draws, seed 20261017). Four of 5,439 nodes retain less than 80% of their members and are untestable. Effects use 367 stage-complete disease clusters (109 GSE162694, 258 GSE213621), adjusted for recorded stage and sex; GSE213621 controls lacking histology are excluded. Systems contain overlapping gene sets and are descriptive summaries.

**(H)** Donor-first fibrosis-stage versus stage-adjusted continuum coefficients for the same 43 systems and 367 clusters, with 95% bootstrap intervals on both axes (1,000 draws). Units are system-score SD per cohort stage-code SD (x) and per continuum SD (y). Colors identify both intervals excluding zero with positive/negative/opposite signs (19/2/1); 21 systems are unresolved. These are nominal interval states, not corrected significance calls. Spearman rho=0.670 (bootstrap 95% interval 0.38–0.84); 34/43 share direction. No maxT or competitive-null outline is displayed. Alignment of these overlapping systems does not establish independent mechanisms or continuum specificity.

See also Figure S4 and Supplementary Notes 8 and 11.

## Figure S4

**Figure S4. Single-cell program, communication and continuum detail, related to Figure 4.**

Bulk panels use selected genetic-identity clusters, not verified independent donors. In continuum displays, dashed gray lines (†) mark GSE126848, GSE130970 and GSE135251, which overlap continuum construction and are descriptive only; solid lines mark GSE162694 and GSE213621. Windows are nine fixed-projection windows, 20 percentile points wide and centered at the 10th to 90th percentiles, so adjacent windows overlap. Window and category means are descriptive. Untestable is not a negative result.

**(A)** Corrected equal-variance sensitivity for Stromal ECM and Ductular injury, with at most 61 source records. Two GSE189600 records have unresolved biological independence. Pooled, cohort and cohort-omission estimates retain their actual denominators and nominal HC3 95% intervals; these are not complete-family primary calls.

**(B)** Corrected MASH-minus-steatosis projected scores in fibroblasts (Stromal ECM) and cholangiocytes (Ductular injury), using the same correction-available cells and donors at the 50-cell gate. Cohort-specific HC3 95% CIs and full donor rosters accompany the source table. The two programs form one BH family per declared arm. Projected scores are not native Hotspot units.

**(C)** Donor-score means after centering each program across the fitted donor universe, shown for Healthy, Steatosis and Steatohepatitis source categories for all 117 fixed programs. The heatmap is descriptive; the accompanying GLS stage beta adjusts for dataset with separate residual variances. Black bars, nominal 95% CI; gray bars, 117-program permutation interval. Source metadata and fitted denominators match Figure 4C.

**(D)** Donor-collapsed communication scores for the fixed 13-pair roster across Healthy, Steatosis and Steatohepatitis (69 donors). Six pairs are testable, all with corrected q = 1. Seven are untestable because stage counts are insufficient; these are drawn at lower opacity and are not tested nulls.

**(E)** Corrected program–communication associations in the 26-test BH family, with 20 estimable and six untestable tests. One Ductular injury association passes both ordinary and HC3 BH: CDH1–PTPRM, 57 donors, beta 0.4984, q 9.940e-6, HC3 q 5.952e-10. Its sign agrees in three of four cohort omissions. Bars are nominal model-based 95% CIs; untestable rows are listed in the complete source table. RNA covariation does not establish signaling.

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

**(V)** Stage- and sex-adjusted and composition-adjusted continuum effects of Stromal ECM and Ductular injury in the 367 disease units (GSE162694, n = 109; GSE213621, n = 258), with composition adjustment under both cell-type fraction sets. Bars are HC3 95% CIs. The weight denominator is retained program weight, not observed weight: Stromal ECM retains 80.5% of its weight but 73.8% is observed. No composition-adjusted program or system is named as a finding.

**(W)** NAS-axis companion to Figure 4E: program-score SD per recorded NAS point for 112 programs, from the same per-cohort OLS fits (sex, fibrosis stage, NAS; HC3) in 364 genetic-identity clusters, combined by REML random-effects meta-analysis with Knapp-Hartung 95% CIs. Colors give discovery lineage. The five largest positive effects are labeled. BH over the 224 tests of Figure 4E; no NAS association meets q < 0.05 (smallest q 0.159).

## Figure 5

**Figure 5. Liver protein tests the disease-state RNA changes, and intact tissue places Figure 5 genes in cell types.**

**(A)** Liver protein against bulk mRNA change for 3,054 measured genes of a 10,044-gene prioritized universe. The universe is partly selected on RNA: rho is 0.487 in RNA-selected genes and 0.305 in the 664 genes selected on genetic evidence alone. Protein model: 46 MASLD versus 12 control participants, adjusted for acquisition batch, age, BMI and sex. Blue, protein-significant (BH p < 0.05 over 6,071 proteins) with the same sign as mRNA; magenta, opposite sign; gray, not significant. The bar shows the 211 protein-significant genes: 188 with the same sign as mRNA and 23 with the opposite sign.

**(B)** 25 proteins fixed before the bulk refit, grouped by process. The same cohort helped select them, so the view is descriptive. Heatmap, protein z-scores, ordered by fibrosis stage then NAS. Middle, partial Spearman rho with each histology score within 46 MASLD participants, participant-bootstrap 95% intervals; filled, BH q < 0.05 over 125 tests. Smaller n are printed beside rows. Right, filled markers are calls (BH padj < 0.05, plus absolute log2FC > 0.5 for mRNA). Process NES from preranked fgsea, BH over 1,145 sets.

**(C)** Cell type that makes each gene in intact healthy liver. Xenium in situ data of Tzouanas et al. (2026; 477-gene panel), 9 control livers; the 9 MASH-cirrhosis explants of the same source are not used here. Rows: the 14 on-panel genes of panel A with a protein change in the same direction as RNA (lower or higher in MASLD), and the 11 on-panel colocalized genes detected in at least 5% of their top cell type's cells in at least 7 of 9 controls (SPI1 is not detected). Dots: spillover-corrected counts per cell, averaged over donors with equal weight, scaled to the gene's top cell type. Endothelial, stellate/smooth-muscle, macrophage and T/NK columns show the highest of their subtypes. Right, the gene's top snRNA lineage in the atlas; "same" when it falls in the top tissue cell type (12 of 14 and 8 of 11). Cell labels and coordinates are the source's. Descriptive; no test.

**(D)** The 22 of 95 scorable fixed programs whose liver protein score is robustly associated with histologic burden in 46 MASLD participants: partial Spearman rho adjusted for acquisition batch, age, BMI and sex, with participant-bootstrap 95% CI. Robust: Freedman-Lane BH q < 0.05 over the 95 scorable programs and stable sign across score variants and both acquisition batches. Rows are selected on the plotted statistic; all 117 programs are in the source table. Blue, lower with greater burden (14); red, higher (8). Histologic burden is the first principal component of steatosis, ballooning, lobular inflammation and fibrosis stage.

## Figure S5: disease accessibility and regulatory-context additions

**Figure S5P–AA. Donor-level disease accessibility, descriptive links and internal predictive controls.** These additions accompany the existing S5 panels. Controls and baseline methods are gray. Accessibility effects are MASH versus normal log2 normalized fragment changes; RNA effects are log2 normalized UMI changes. Neither peak significance nor predictive performance validates a linked gene.

**(P)** ATAC nuclei per donor in the two disease cohorts and three primary lineages. Normal, gray; MASH, magenta. Donors, rather than nuclei, are the inferential units. The primary accessibility and the 100-nucleus network eligibility thresholds are distinct.

**(Q,R)** Cohort accessibility effects for all 30,650 jointly testable hepatocyte peaks and 9,264 stellate peaks, respectively. Spearman rho 0.390 and 0.265 are descriptive across peaks, which are not independent donors. Four stellate peakwise candidates are highlighted; no lineage passes the separate three-lineage directional-replication BH family.

**(S)** Condition-only versus recorded-batch accessibility effects in GSE244832 on identical 14 donors and fixed peaks. Counts printed on the panel are BH q < 0.05 within lineage and model: hepatocytes 1 versus 3, stellate 820 versus 291. Gray diagonal denotes equal estimates. Batch and disease are not fully separable merely because the design is full rank.

**(T)** Each donor's omission effect on the median absolute change in accessibility log2 effect, by cohort and lineage. These are independently refitted donor omissions, not confidence intervals.

**(U)** Cross-cohort Jaccard of selected peak–gene links among pairs observable in both cohorts: Cicero 0.158 and 0.174, distance 0.993 and 0.998, for hepatocytes and stellate cells. Shared consensus ascertainment and fixed geometry qualify this comparison; stability is not regulatory accuracy.

**(V)** Original L1 program weight with a positive Cicero annotation, retaining all 117 fixed programs in each cohort/lineage context. Link coverage is distinct from program significance; memberships and weights are unchanged.

**(W)** Deterministically selected accessibility examples, including a peakwise candidate and source-dependent, indeterminate and batch-sensitive contexts. Points and bars are peak effects and pointwise 95% moderated-t intervals. Gene names indicate descriptive links only. The example-selection rule used the complete evidence table, not expected mechanisms.

**(X)** Cicero RNA-prediction gain over retuned promoter-only and distance models in 39 GSE296875 held-out donors. Positive is lower standardized MSE. Points and bars show donor-average differences and 95% intervals from 2,000 seeded donor bootstraps conditional on fitted models. The original five outer/three inner donor folds and wider ridge grid are fixed. Aggregation seeds reuse the original capped nuclei; omissions are targeted training-only stress tests. These descriptive intervals add no confirmatory family and omit retraining uncertainty.

**(Y)** Correctly paired minus RNA-shuffled Cicero gain over distance for three within-donor/lineage shuffles. All conditional 95% donor intervals include zero. Shuffles preserve donor RNA marginals and do not remove all donor-level association.

**(Z)** Hepatocyte RNA candidates CHST9 and TMEM150B, which meet BH q < 0.05 in both condition-only and recorded-batch models, and TRIM9, which meets the threshold only after adjustment. The displayed genes are selected from the tested family; pointwise 95% moderated-t intervals are not simultaneous or selection-adjusted. Filled markers meet gene-level BH within model across 21,054 genes; open markers do not. Both models use 5 normal and 9 MASH donors. The complete family is in the RNA effects table. These are separate-library, single-cohort associations, with no new ambient correction or doublet model.

**(AA)** Probability of any peak-level BH rejection in 1,000 fitted-negative-binomial null simulations for fragment counts, with Monte Carlo uncertainty. Gray, primary voom; magenta, robust edgeR QL. This model-null diagnostic is not biological FDR. Binary accessible-nucleus simulations are excluded because negative-binomial counts do not reproduce their bounded support.

## Figure S5: panels moved from the main figure (2026-11-03; letters to be assigned)

These were Figure 5A, 5D, 5E and 5F before 2026-11-03. Their legends are kept unchanged below, with the former main-figure letter. In the former 5E legend, "panel D" means the former 5D.

**(former 5A)** Inputs and readouts. Each input set is fixed before its readout; no program is reselected.

**(former 5D)** Percent of colocalized genes whose best gene-study pair places at least 50% of its primary-signal posterior in open chromatin of each lineage, in at least one of two liver snATAC sources (circle) or both (triangle). Direct MASLD/liver-fat genes, n = 22 (28 pairs); liver-enzyme genes, n = 424 (739 pairs); ten genes are in both. Intervals, gene-bootstrap 95%. Diamonds, expectation with posterior spread uniformly over tested variants; all ten are 0% by construction, because an evenly spread posterior cannot place 50% in peaks. Against a matched expectation, the observed values are 1.8 to 2.4 times expected, not 2.5 to 3.1. Only hepatocyte-specific peaks are enriched (1.77); stellate, macrophage and cholangiocyte-specific peaks are below expectation (0.42 to 0.62). Overlap does not assign a target gene. Disease-associated accessibility is evaluated separately in Supplementary Figure S5Q–T.

**(former 5E)** 45 of 117 programs; all 117 are in the source table. Upper block, 31 programs with a member gene open in its own lineage through a panel D pair. Lower block, 14 programs shown only for a robust protein association, so selected on the plotted statistic. Bars: dark, open through a direct MASLD/liver-fat pair (RORA); light, through liver-enzyme pairs only. Right, partial Spearman rho (participant-bootstrap 95% CI) between program protein score and histologic burden in 46 MASLD participants. Robust: Freedman-Lane BH q < 0.05 over 95 scorable programs and stable sign across score variants and both batches. Gray outline, CI excludes 0 but not robust.

**(former 5F)** Left, illustrative H&E overlay of the adjusted Stromal ECM score (section JBO019, donor H37), chosen outcome-blind; no pathologist annotation. Right, Moran excess over matched genes after adjusting spot scores for all 16 cell2location abundances, in Guilliams et al. and Vu et al. (2026). Small points, donors (Guilliams et al., 4) or arrays (Vu et al., 10); band, matched-null 95% range. Labels give BH q from 9,999 matched-gene draws (two programs per source) and donor or array counts. The gene null holds donors and arrays fixed and draws control genes independently, matched on correlation with hepatocyte abundance, so it is not a donor-level test and is not matched for program co-expression. In a descriptive analysis that was not prespecified, 68 of 114 (Guilliams et al.) and 74 of 114 (Vu et al.) scorable fixed programs exceed a null built the same way, with each program's controls matched on its own discovery lineage (BH q < 0.05 over 114 programs per source); Stromal ECM ranks 42nd and 29th by excess. With controls matched on correlation with fibroblast abundance, Stromal ECM still exceeds the null (q = 0.011 and 0.0008); with controls matched on correlation with cholangiocyte abundance, Ductular injury does not (q = 0.10 and 0.50). Donor and array counts are of units with positive excess over the matched null: Stromal ECM excess is positive in 3 of 4 Guilliams et al. donors (exact one-sided sign-flip p = 0.19; the smallest possible p with 4 donors is 0.0625) and Ductular injury excess in 3 of 4 (p = 0.125). Spot scores use the measured member genes, which carry 55% (Guilliams et al.) and 53% (Vu et al.) of the Stromal ECM L1 weight and 41% of the Ductular injury weight. IGFBP7 carries 2.47% of the full Stromal ECM weight and BICC1 5.44% of the full Ductular injury weight, so the result concerns program scores, not either gene. The two sources used different cell2location reference builds (Methods). Key: "supported" means the program exceeds the matched-gene null. Vu et al. has no donor key (open point, source-dependent; this applies to Stromal ECM only). The Guilliams et al. donors are surgical samples taken for other conditions, not a MASLD cohort; H37 and H38 are also in the single-cell atlas.

## Figure 6

**Figure 6. The MASLD Gene Catalog keeps eQTL-colocalization and bulk RNA calls separate for 28,233 genes.**

**(A)** Catalog structure. Each gene identity carries a separate eQTL-colocalization call and bulk RNA call, each qualified by testability and provenance, plus one molecular object (protein-coding gene, lncRNA transcript, ambiguous lncRNA locus or untestable nomination; no record is a DNA element) and one proposed discriminating experiment. Protein, chromatin and spatial context keep their assay or program unit and enter no class or route.

**(B)** Genetic call (rows) by bulk RNA call (columns), 28,233 identities, in the November 6 eligibility release. Genetic support requires SuSiE-COLOC PP.H4 > 0.5 in any of 35 Tier-1/2 GWAS; RNA support requires BH q < 0.05 over 24,115 genes and |log2FC| > 0.5 in 774 genetic-identity clusters (138 control, 636 disease), not verified donors. The genetic-state rows supported, indeterminate and untestable contain 26/400/10, 1,118/12,958/4,214 and 709/8,790/8 genes across RNA columns supported/indeterminate/untestable. Blue means neither domain is classed untestable; gray is untestable in at least one. Indeterminate genetics can include an ABF record without a SuSiE result and is not a negative multi-signal test. Genetic untestability includes 9,312 identities without a primary COLOC record and 195 with a priority untestable state, for 9,507 total. RNA is untestable for 4,232; eight are untestable in both. No negative-call rule exists, so no tested-negative row or column is shown. Coding PIP ≥ 0.5 annotates 31 genes and does not enter the eQTL-colocalization call.

**(C)** Four worked entries. Columns give SuSiE and ABF PP.H4 with the top study, any fine-mapped coding variant, bulk RNA log2FC with 95% CI (log2FC ± 1.96 SE; 774 genetic-identity clusters) with BH q, and the Catalog class. The line under each gene gives its gene-level and molecular-object routes. Dashed lines mark |log2FC| = 0.5. ABF PP.H4 is shown for reference and makes no call. *GLP1R* is untested because its RNA was not detected; its genetic call is indeterminate (no Tier-1/2 SuSiE result). *MTARC1* has no Tier-1/2 SuSiE result, so its genetic call is indeterminate. *MTARC1* and *THRB* have indeterminate RNA calls because the change is below the 0.5 floor. The *MTARC1* missense variant (PIP 0.76) does not enter its genetic call. "Both supported" (*HKDC1*) does not compare direction. Genetic support does not establish a causal mechanism.

**(D)** Genes per adopted Catalog class, with the gene-level route: untested 13,731, unresolved 12,958, established-state associated 1,118, genetically anchored 400 and supported on both axes 26 (class concordant, without comparing directions). The untested bar partitions 4,232 RNA-untestable genes, including eight untestable in both domains, and 9,499 genetics-untestable genes whose RNA is testable. Gray is the untested class; blue shows the other classes. Hatching marks 719 untested genes with one supported call (709 RNA, 10 genetics). Twenty unresolved genes carry a coding variant at PIP ≥ 0.5 but keep their eQTL-colocalization class. Among unresolved, established-state and untested genes, 290, 18 and 25 have ABF PP.H4 > 0.5 without a primary SuSiE result. The November 6 eligibility qualification changes 180 classes without changing positive genetic genes or RNA calls. Routes are proposals, not a priority ranking.

**(E)** Ways to use the Catalog: view the worked entries, inspect each molecular object with its uncertainty, read source references, and download the gene calls with the molecular-object records. All views draw on the same inputs and give no combined score or rank. The Catalog files are available at [accession to be supplied] (see Data availability).

## Figure 7

**Figure 7. Predictive models are scored against task-specific baselines beside a record of their training exposure.**

**(A)** Evaluation framework. Each of 8 registered tasks has a baseline set in advance. A score whose exposure record is unresolved is uninterpretable, not evidence of skill. The bar counts 126 audit entries by checkpoint finding, including 37 no-finding placeholders, before eight later audit revisions. These checkpoint findings are separate from the study-specific exposure in (C) and (F).

**(B)** Gain over each task's own baseline, named with its score on each x axis. Reporter assay (GSE281364, HepG2 contexts only): Bonferroni-corrected intervals across the seven rows (paired block bootstrap). Sequence-based reporter baselines were run: HyenaDNA scores 0.216 against allele-identity ridge (B0) 0.134, CpG ridge (B1) 0.178, composition (B2) 0.176 and the prespecified flank and dinucleotide ridge (B3) 0.167. The gain over B3 is 0.049 (95% 0.0006 to 0.099; Bonferroni-7 lower bound -0.015) and over B1 0.038 (95% -0.006 to 0.083), so the corrected gain holds only against allele identity. eQTL |z| and sign: held-out LD fold, paired 95% LD-block bootstrap. Fibrosis: samples are the unit because the training pool has no donor key; 95% sample bootstrap. H3K27ac: 95% participant bootstrap. ATAC: donor jackknife.

**(C)** Donor-class-balanced macro F1 for broad cell classes on all seven held-out studies (hollow; 102 donors, 50,000 cells) and on each encoder's clean studies (filled). Gray: HVG+PCA with the same head. Intervals are 95% study-stratified donor-cluster bootstraps. Geneformer V2 has no clean study, so its scores are not evidence of skill. Hollow and filled scores use different study sets and are not compared with each other; (F) compares each encoder with HVG+PCA on the same studies.

**(D)** Signed eQTL: best head's sign AUC per |z| stratum (95% LD-block bootstrap, uncorrected across five strata) beside within-block permutation chance and the sequence-attributable bound, which assumes at most one causal variant per gene and LD block. Severity transfer: released instrument, no refit; Spearman with stage (stage-stratified 95% bootstrap) beside each cohort's stage tie ceiling. GSE213621 units are genetic identities, not verified donors. GSE276114 includes participants with HCC, sampled from peritumoral liver without a per-sample label. LD-block join: Pearson r of each row's effect with its gene's mean effect, own LD block left in or held out. Routing: fitted rule and oracle router against the pre-declared floor.

**(E)** Each statistic is oriented so that positive is the prespecified direction and placed on one scale, (estimate − null)/(threshold − null): 0 is the null, 1 the threshold its verdict uses (minimum detectable effect for severity arms and the zone task; h2 z = 7 for GWAS rows; twice the test's resolution for signed eQTL; the pre-declared floor for routing). Severity-arm intervals (95% stratified bootstrap) show precision only; verdicts use the label permutation null, and an underpowered arm is untestable, not negative. Signed eQTL direction, recorded as tested negative, is drawn untestable because its sequence-attributable bound (open diamond) is below the threshold. Treatment response uses the GSE83452 diet arm (6 resolved vs 16 persisted; permutation p 0.23). Its underpowered result is drawn as untestable rather than a tested negative. Routing is also untestable because its comparison does not establish an adequately powered negative result.

**(F)** Top: exposure of each held-out study (donors in parentheses) for each encoder. Middle: macro F1 gain over HVG+PCA with the same head, on all seven studies (hollow) and clean studies (filled). Bottom: encoder minus the scVI liver encoder trained here, on the encoder's clean studies, drawn as the donor-pooled estimate. For TranscriptFormer minus scVI (linear head), the donor-pooled difference is +0.0154 and the equal-study mean is -0.0042; both intervals exclude zero, so the sign depends on the estimand. Intervals are 95% donor-cluster bootstraps, uncorrected across rows. In the unweighted study mean, a study with 2 donors and 149 cells gets one seventh of the weight; two of twelve model-by-head blocks change sign when weighted by donors, so no sign claim rests on the unweighted mean.
