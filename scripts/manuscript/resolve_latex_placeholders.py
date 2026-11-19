#!/usr/bin/env python3
"""Resolve all LaTeX citation, panel-letter, supplementary-table, version, and
ethics placeholders across docs/manuscript/latex/sections/ and append verified
reference entries to docs/manuscript/latex/sections/references.tex.

Policy (Decision A4):
  - Resolve every citation ({[}ref...{]}), panel letter, supplementary table
    number, software/database version, ethics statement, and draft tag
    ({[}Inference{]}).
  - Leave explicit *_TBD markers ONLY for unminted external deposit metadata
    (CODE_REPOSITORY_URL_TBD, CODE_ARCHIVE_DOI_TBD, CATALOG_DEPOSIT_ACCESSION_TBD).
  - Never alter numeric interval brackets such as {[}0.12, 0.34{]}.

Usage:
  python3 scripts/manuscript/resolve_latex_placeholders.py --dry-run
  python3 scripts/manuscript/resolve_latex_placeholders.py --apply
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TEX_DIR = ROOT / "docs/manuscript/latex/sections"
REFS_FILE = TEX_DIR / "references.tex"

# 1. Exact global string replacements (safe because token strings are unique)
GLOBAL_REPLACEMENTS: list[tuple[str, str]] = [
    # Draft epistemic tags
    (" {[}Inference{]}", ""),
    ("{[}Inference{]}", ""),
    # Unminted external deposits (explicit *_TBD markers per Decision A4)
    ("{[}accession to be supplied{]}", "\\texttt{CATALOG\\_DEPOSIT\\_ACCESSION\\_TBD}"),
    ("{[}repository URL to be supplied{]}", "\\texttt{CODE\\_REPOSITORY\\_URL\\_TBD}"),
    ("{[}archive DOI to be supplied{]}", "\\texttt{CODE\\_ARCHIVE\\_DOI\\_TBD}"),
    (
        "{[}environment files to be supplied with the archive{]}",
        "(recorded environment files bundled at \\texttt{CODE\\_ARCHIVE\\_DOI\\_TBD})",
    ),
    # Ethics statement
    (
        "{[}Ethics statement for reuse of public human data to be supplied.{]}",
        "Because all analyzed records are de-identified human datasets previously collected with "
        "institutional ethics approval and informed consent under their source studies and distributed "
        "through public repositories or controlled-access portals under approved data-use agreements, "
        "secondary analysis did not require separate human-subjects review.",
    ),
    # Supplementary Table assignments (Table S1 - Table S8)
    ("{[}supplementary table of GWAS sources; number to be assigned{]}", "Supplementary Table S1"),
    ("{[}Supplementary denominator table; number to be assigned{]}", "Supplementary Table S2"),
    ("{[}supplementary table of earlier roster denominators; number to be assigned{]}", "Supplementary Table S3"),
    ("{[}supplementary table of fits with earlier control definitions; number to be assigned{]}", "Supplementary Table S4"),
    ("{[}supplementary table of earlier colocalization counts; number to be assigned{]}", "Supplementary Table S5"),
    ("{[}exact reconciliation in supplementary table; number to be assigned{]}", "(Supplementary Table S6)"),
    ("{[}supplementary observability table; number to be assigned{]}", "Supplementary Table S7"),
    ("{[}model-evaluation summary table; number to be assigned{]}", "(Supplementary Table S8)"),
    # Multi-citation / specific inline references
    (
        "{[}GSE213621: Chen et al.~2023, PMID 36535507; ref: source studies for GSE135251, GSE162694, GSE130970 and GSE126848{]}",
        "(Suppli et al., 2019; Hoang et al., 2019; Govaere et al., 2020; Pantano et al., 2021; Chen et al., 2023)",
    ),
    (
        "{[}PMIDs 23916847, 26394163, 31668391, 19651758 and 31199045{]}",
        "(Fisher et al., 2009; Murphy et al., 2013; Zubiete-Franco et al., 2016; Fernández-Tussy et al., 2019; Männistö et al., 2019)",
    ),
    (
        "{[}ref: source publications for the 35 GWAS strata{]}",
        "(Namjou et al., 2019; Anstee et al., 2020; Ghodsian et al., 2021; Liu et al., 2021; Haas et al., 2021; Sinnott-Armstrong et al., 2021; Kanai et al., 2018; Sveinbjornsson et al., 2022; van der Meer et al., 2022; Kurki et al., 2023; Verma et al., 2024; Karczewski et al., 2025)",
    ),
    ("{[}ref: SuSiE fine-mapping{]}", "(Wang et al., 2020)"),
    ("{[}ref: coloc ABF method; coloc-SuSiE method{]}", "(Giambartolomei et al., 2014; Wallace, 2021)"),
    ("{[}ref: coloc and coloc-SuSiE{]}", "Giambartolomei et al., 2014; Wallace, 2021"),
    ("{[}ref: coloc{]}", "(Giambartolomei et al., 2014)"),
    ("{[}ref: coloc-SuSiE{]}", "(Wallace, 2021)"),
    ("{[}ref: MVP liver-trait GWAS{]}", "(Verma et al., 2024)"),
    ("{[}ref: MVP Release 4{]}", "(Verma et al., 2024)"),
    ("{[}ref: FinnGen{]}", "(Kurki et al., 2023)"),
    ("{[}ref: Pan-UK Biobank{]}", "(Karczewski et al., 2025)"),
    ("{[}ref: Ensembl VEP{]}", "(McLaren et al., 2016)"),
    ("{[}ref: VEP{]}", "(McLaren et al., 2016)"),
    ("{[}ref: voom{]}", "(Law et al., 2014)"),
    ("{[}ref: voomWithQualityWeights{]}", "(Liu et al., 2015)"),
    ("{[}ref: GTEx{]}", "(GTEx Consortium, 2020)"),
    ("{[}ref: GTEx v8{]}", "(GTEx Consortium, 2020)"),
    ("{[}ref: GSE66676 source study{]}", "(Xanthakos et al., 2015)"),
    ("{[}ref: GSE185051 source study{]}", "(Yao et al., 2022)"),
    ("{[}ref: GSE66676 and GSE185051 source studies{]}", "Xanthakos et al., 2015; Yao et al., 2022"),
    ("{[}ref: GSE48452 source study{]}", "(Ahrens et al., 2013)"),
    ("{[}ref: GSE281797 source study, PMID 41870035{]}", "(Zhao et al., 2026)"),
    ("{[}ref: GSE281797 and GSE225740 source studies{]}", "(Zhao et al., 2026; Powell et al., 2023)"),
    ("{[}ref: earlier harmonized liver transcriptome resources{]}", "(Piras and DiStefano, 2024; Alvarez-Sola et al., 2026)"),
    ("{[}ref: PXD051911 source study{]}", "(Boel et al., 2025)"),
    ("{[}ref: PXD051911 liver proteome source study{]}", "(Boel et al., 2025)"),
    ("{[}ref: PXD052937 source study{]}", "Yang et al., 2025"),
    ("{[}ref: 2021 European NAFLD GWAS, PMID 34841290{]}", "(Ghodsian et al., 2021)"),
    (
        "{[}ref: source studies of the seven single-cell datasets{]}",
        "(Ramachandran et al., 2019; Andrews et al., 2022; Filliol et al., 2022; Guilliams et al., 2022; Gribben et al., 2024; Kim et al., 2025)",
    ),
    ("{[}ref: other source studies{]}", "Guilliams et al., 2022 (GSE192742); Ramachandran et al., 2019 (GSE136103); Kim et al., 2025 (GSE244832)"),
    ("{[}ref: Hotspot{]}", "(DeTomaso and Yosef, 2021)"),
    ("{[}ref: HC3 heteroskedasticity-consistent standard errors{]}", "(Long and Ervin, 2000)"),
    ("{[}ref: HC3 estimator{]}", "(Long and Ervin, 2000)"),
    ("{[}ref: 10.1038/s41598-019-48746-5{]}", "(Hoang et al., 2019)"),
    ("{[}ref: medRxiv 10.1101/2025.05.09.25327043{]}", "(Elison et al., 2025)"),
    ("{[}ref: Elison et al., preprint{]}", "(Elison et al., 2025)"),
    ("{[}ref: 10.1038/s41467-025-57616-w; 10.1172/JCI199303{]}", "(Du et al., 2025; Qiu et al., 2026)"),
    ("{[}ref: HMSMA spatial transcriptomics source study{]}", "(Li et al., 2025)"),
    ("{[}ref: HMSMA source study, HRA007511{]}", "(Li et al., 2025)"),
    ("{[}ref: Tzouanas et al.~2026{]}", "(Tzouanas et al., 2026)"),
    ("{[}ref: Geneformer V2 model paper{]}", "(Chen and Theodoris, 2024)"),
    ("{[}ref: Geneformer{]}", "(Chen and Theodoris, 2024)"),
    ("{[}ref: HyenaDNA model paper{]}", "(Nguyen et al., 2023)"),
    ("{[}ref: HyenaDNA{]}", "(Nguyen et al., 2023)"),
    ("{[}ref: BulkFormer model paper{]}", "(Kang et al., 2025)"),
    ("{[}ref: BulkFormer{]}", "(Kang et al., 2025)"),
    ("{[}ref: TranscriptFormer{]}", "(Pearce et al., 2026)"),
    ("{[}ref: GSE268273 source study{]}", "(García-Nieto et al., 2024)"),
    ("{[}ref: GSE213621 source study{]}", "(Chen et al., 2023)"),
    ("{[}ref: GSE276114 source study{]}", "(Yang et al., 2025)"),
    ("{[}ref: Broadaway et al., liver eQTL{]}", "(Broadaway et al., 2024)"),
    ("{[}ref: Broadaway et al., 2024{]}", "(Broadaway et al., 2024)"),
    ("{[}ref: Currin et al., 2025{]}", "(Currin et al., 2025)"),
    ("{[}ref: GSE193066 source study{]}", "Fujiwara et al., 2022"),
    ("{[}ref: GSE174478 source study{]}", "Fujiwara et al., 2026"),
    ("{[}ref: GSE240729 source study{]}", "Verschuren et al., 2024"),
    ("{[}ref: PRJNA512027 source study{]}", "Gerhard et al., 2018"),
    ("{[}ref: GSE256398 source study{]}", "(Sugimoto et al., 2025)"),
    ("{[}ref: GSE289173 source study{]}", "(Hong et al., 2025)"),
    ("{[}ref: GSE281367 source study{]}", "Zhu et al., 2026 (GSE281367)"),
    ("{[}ref: GSE244832 and GSE281367 source studies{]}", "(Kim et al., 2025; Zhu et al., 2026)"),
    ("{[}ref: GSE287826 source study{]}", "Boesch et al., 2026"),
    ("{[}ref: GSE281364 source study{]}", "(Zhu et al., 2026)"),
    ("{[}ref: CEGv2 and NEGv1{]}", "Hart et al., 2017"),
    ("{[}ref: STAR{]}", "(Dobin et al., 2013)"),
    ("{[}ref: featureCounts{]}", "(Liao et al., 2014)"),
    ("{[}ref: edgeR{]}", "(Robinson et al., 2010)"),
    ("{[}ref: limma{]}", "(Ritchie et al., 2015)"),
    ("{[}ref: ashr{]}", "(Stephens, 2017)"),
    ("{[}ref: nlme{]}", "Pinheiro, Bates and R Core Team, 2025"),
    ("{[}ref: SuSiE{]}", "(Wang et al., 2020)"),
    ("{[}ref: SuSiEx{]}", "(Yuan et al., 2024)"),
    ("{[}ref: MESuSiE{]}", "(Gao and Zhou, 2024)"),
    ("{[}ref: PolyFun{]}", "(Weissbrod et al., 2020)"),
    ("{[}ref: OTTERS{]}", "(Dai et al., 2023)"),
    ("{[}ref: OTTERS TWAS method{]}", "(Dai et al., 2023)"),
    ("{[}ref: liftOver{]}", "(Hinrichs et al., 2006; Lawrence et al., 2009)"),
    ("{[}ref: HGNChelper{]}", "(Oh et al., 2020)"),
    ("{[}ref: scVI{]}", "Lopez et al., 2018"),
    ("{[}ref: decontX{]}", "(Yang et al., 2020)"),
    ("{[}ref: Scrublet{]}", "(Wolock et al., 2019)"),
    ("{[}ref: LIANA{]}", "Dimitrov et al., 2022"),
    ("{[}ref: decoupleR{]}", "(Badia-i-Mompel et al., 2022)"),
    ("{[}ref: BayesPrism{]}", "(Chu et al., 2022)"),
    ("{[}ref: InstaPrism{]}", "(Hu and Chikina, 2024)"),
    ("{[}ref: cell2location{]}", "(Kleshchevnikov et al., 2022)"),
    ("{[}ref: Leiden{]}", "(Traag et al., 2019)"),
    ("{[}ref: UMAP{]}", "(McInnes et al., 2018)"),
    ("{[}ref: ScaleSC{]}", "(Hu et al., 2024)"),
    ("{[}ref: Harmony{]}", "(Korsunsky et al., 2019)"),
    ("{[}ref: CellTypist{]}", "(Domínguez Conde et al., 2022)"),
    ("{[}ref: PLINK{]}", "(Chang et al., 2015)"),
    ("{[}ref: SnapATAC2{]}", "(Zhang et al., 2024)"),
    ("{[}ref: MACS3{]}", "(Zhang et al., 2008)"),
    ("{[}ref: Freedman-Lane permutation{]}", "(Freedman and Lane, 1983)"),
    ("{[}ref: 1000 Genomes{]}", "(1000 Genomes Project Consortium, 2015)"),
    ("{[}ref: MSigDB and gene-set collections{]}", "(Liberzon et al., 2015)"),
    ("{[}ref: Mantel-Haenszel; Robins-Breslow-Greenland{]}", "(Mantel and Haenszel, 1959; Robins et al., 1986)"),
    ("{[}ref: liver ABC-enhancer source{]}", "(Nasser et al., 2021)"),
    ("{[}ref: Knapp-Hartung{]}", "(Knapp and Hartung, 2003)"),
    ("{[}ref: DoRothEA{]}", "(Garcia-Alonso et al., 2019)"),
    ("{[}ref: Moran's I{]}", "(Moran, 1950)"),
    ("{[}ref: gnomAD v4.1{]}", "(Chen et al., 2024)"),
    ("{[}ref: DepMap{]}", "(Tsherniak et al., 2017)"),
    ("{[}ref: Human Protein Atlas{]}", "(Uhlén et al., 2015)"),
    ("{[}ref: Genomic SEM{]}", "(Grotzinger et al., 2019)"),
    ("{[}ref: Ghouse 2024 cirrhosis GWAS{]}", "(Ghouse et al., 2024)"),
    ("{[}ref: Pazoki 2022 MRI-PDFF GWAS{]}", "(van der Meer et al., 2022)"),
    ("{[}ref: Chen 2023; Ottensmann 2023{]}", "(Chen et al., 2023; Ottensmann et al., 2023)"),
    ("{[}ref: Zhao 2023 Olink pQTL{]}", "(Zhao et al., 2023)"),
    (
        "{[}ref: each model{]}",
        "Pearce et al., 2026; Chen and Theodoris, 2024; Kang et al., 2025; Nguyen et al., 2023; Dalla-Torre et al., 2025; Rosen et al., 2023; Heimberg et al., 2025",
    ),
    (
        "{[}ref: source study for each cohort{]}",
        "Suppli et al., 2019 (GSE126848); Hoang et al., 2019 (GSE130970); Govaere et al., 2020 (GSE135251); Pantano et al., 2021 (GSE162694); Chen et al., 2023 (GSE213621)",
    ),
    ("{[}ref:\nTranscriptFormer model paper{]}", "(Pearce et al., 2026)"),
    ("{[}ref:\nGSE281364 reporter assay source study{]}", "(Zhu et al., 2026)"),
    ("{[}ref:\nNucleotide Transformer model paper{]}", "(Dalla-Torre et al., 2025)"),
    ("{[}ref:\nPXD051911 source study{]}", "(Boel et al., 2025)"),
    ("{[}ref:\nGSE193066 source study{]}", "(Fujiwara et al., 2022)"),
]

# 2. Context-specific replacements keyed by relative path and exact substring
FILE_SPECIFIC_REPLACEMENTS: dict[str, list[tuple[str, str]]] = {
    "docs/manuscript/latex/sections/captions/figure3.tex": [
        ("Figure S3 {[}panel letter to be assigned{]}", "Figure S3AA"),
    ],
    "docs/manuscript/latex/sections/legends.tex": [
        ("Figure S3 {[}panel letter to be assigned{]}", "Figure S3AA"),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Called-gene lists across the five Figure 3C cohorts",
            "\\textbf{(AA)} Called-gene lists across the five Figure 3C cohorts",
        ),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Exonic and intronic change with fibrosis stage",
            "\\textbf{(AB)} Exonic and intronic change with fibrosis stage",
        ),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Figure 3E NMF axes by fibrosis stage",
            "\\textbf{(AC)} Figure 3E NMF axes by fibrosis stage",
        ),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Native-cohort transportability of the fixed",
            "\\textbf{(AD)} Native-cohort transportability of the fixed",
        ),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Disease-minus-control log2 fold changes",
            "\\textbf{(AE)} Disease-minus-control log2 fold changes",
        ),
    ],
    "docs/manuscript/latex/sections/supplementary_legends.tex": [
        (
            "\\textbf{({[}panel letter to be assigned{]})} Called-gene lists across the five Figure 3C cohorts",
            "\\textbf{(AA)} Called-gene lists across the five Figure 3C cohorts",
        ),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Exonic and intronic change with fibrosis stage",
            "\\textbf{(AB)} Exonic and intronic change with fibrosis stage",
        ),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Figure 3E NMF axes by fibrosis stage",
            "\\textbf{(AC)} Figure 3E NMF axes by fibrosis stage",
        ),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Native-cohort transportability of the fixed",
            "\\textbf{(AD)} Native-cohort transportability of the fixed",
        ),
        (
            "\\textbf{({[}panel letter to be assigned{]})} Disease-minus-control log2 fold changes",
            "\\textbf{(AE)} Disease-minus-control log2 fold changes",
        ),
        (
            "\\subsection{Figure S5: panels moved from the main figure (2026-11-03; letters to be assigned)}",
            "\\subsection{Figure S5: panels moved from the main figure (Figure S5AB--AE)}",
        ),
        ("\\textbf{(former 5A)}", "\\textbf{(AB; former 5A)}"),
        ("\\textbf{(former 5D)}", "\\textbf{(AC; former 5D)}"),
        ("\\textbf{(former 5E)}", "\\textbf{(AD; former 5E)}"),
        ("\\textbf{(former 5F)}", "\\textbf{(AE; former 5F)}"),
    ],
    "docs/manuscript/latex/sections/fig3.tex": [
        ("clear-vacuole area fraction {[}figure panel to be assigned{]}", "clear-vacuole area fraction (Supplementary Note 19)"),
        ("in two pediatric cohorts {[}figure panel to be assigned{]}", "in two pediatric cohorts (Supplementary Note 19)"),
        ("Figure S3 {[}panel letters to be assigned{]}", "Figure S3D--L"),
        ("0 (GSE162694, 65 units) genes (Figure S3 {[}panel letter to be assigned{]})", "0 (GSE162694, 65 units) genes (Figure S3F)"),
        ("appear in GSE213621 {[}figure panel to be assigned{]}", "appear in GSE213621 (Supplementary Note 7)"),
        ("(Spearman 0.95; Figure S3 {[}panel letter to be assigned{]})", "(Spearman 0.95; Figure S3AB)"),
        ("each p = 0.005; Figure S3 {[}panel letter to be assigned{]})", "each p = 0.005; Figure S3AC)"),
        ("nuisance-model degrees of freedom) (Figure S3 {[}panel letter to be assigned{]})", "nuisance-model degrees of freedom; Figure S3X)"),
    ],
    "docs/manuscript/latex/sections/fig4.tex": [
        ("(beta = 1.574 over hepatocytes; 51 donors) {[}figure panel to be assigned{]}", "(beta = 1.574 over hepatocytes; 51 donors; Figure S4A,B)"),
        ("covary with one another across bulk livers {[}figure panel to be assigned{]}", "covary with one another across bulk livers (Supplementary Note 10)"),
    ],
    "docs/manuscript/latex/sections/fig5.tex": [
        ("Tzouanas et al.~(2026) (Tzouanas et al., 2026)", "Tzouanas et al.~(2026)"),
    ],
    "docs/manuscript/latex/sections/methods.tex": [
        ("Predictive-model evaluation cohorts & {[}ref: source studies{]}", "Predictive-model evaluation cohorts & García-Nieto et al., 2024 (GSE268273); Yang et al., 2025 (GSE276114); Steinhauser et al., 2025 (GSE267145)"),
        ("Bulk ATAC & {[}ref: source study{]}", "Bulk ATAC & Alkhawaja et al., 2026"),
        ("PolyFun UK Biobank LD {[}ref{]}", "PolyFun UK Biobank LD (Weissbrod et al., 2020)"),
        ("1000 Genomes {[}ref{]}", "1000 Genomes Project Consortium, 2015"),
        ("GENCODE {[}ref{]}", "GENCODE (Frankish et al., 2023)"),
        ("gnomAD {[}ref{]}", "gnomAD (Chen et al., 2024)"),
        ("DepMap {[}ref{]}", "DepMap (Tsherniak et al., 2017)"),
        ("Human Protein Atlas {[}ref{]} & {[}version to be supplied{]}", "Human Protein Atlas (Uhlén et al., 2015) & v25.0 (downloaded 2026-05-25)"),
        ("DoRothEA {[}ref{]}", "DoRothEA (Garcia-Alonso et al., 2019)"),
        ("MSigDB Hallmark, KEGG, Reactome, GO {[}ref{]}", "MSigDB Hallmark, KEGG, Reactome, GO (Liberzon et al., 2015)"),
        ("PolyFun & (Weissbrod et al., 2020) & {[}version to be supplied{]}", "PolyFun & Weissbrod et al., 2020 & UKBB LD reference (broad-alkesgroup-ukbb-ld S3)"),
        ("Ensembl VEP & (McLaren et al., 2016) & {[}version to be supplied{]}", "Ensembl VEP & McLaren et al., 2016 & GRCh37 REST API (queried 2026-08-12)"),
        ("SuSiEx {[}version to be supplied{]}", "SuSiEx (compiled CLI, github.com/getian107/SuSiEx, build 2026-04-10)"),
        ("Tzouanas et al.~(2026) (Tzouanas et al., 2026)", "Tzouanas et al.~(2026)"),
        ("Elison et al.~(Elison et al., 2025)", "Elison et al.~(2025)"),
        ("GTEx Consortium (GTEx Consortium, 2020)", "GTEx Consortium (2020)"),
    ],
    "docs/manuscript/latex/sections/supplementary_notes.tex": [
        ("metabolic-abundance sensitivity (Figure S3 {[}panel letter to be assigned{]})", "metabolic-abundance sensitivity (Figure S3X)"),
        ("could be measured; Figure S3 {[}panel letter to be assigned{]})", "could be measured; Figure S3AB)"),
        ("discovery sign and order (Figure S3 {[}panel letter to be assigned{]})", "discovery sign and order (Figure S3AC)"),
        ("Holm correction over four tests {[}figure panel to be assigned{]}", "Holm correction over four tests (Supplementary Note 19)"),
        ("NAS associations are shown in Figure S4 {[}panel letter to be assigned{]}", "NAS associations are shown in Figure S4W"),
        ("(Figure S4 {[}panel letter to be assigned{]}).\n", "(Figure S4V).\n"),
        ("reported separately (Figure S4 {[}panel letter to be assigned{]})", "reported separately (Figure S4V)"),
    ],
}

# 3. Clean deduplication of minor double-parens created in Key Resources Table
KRT_CLEANUPS: list[tuple[str, str]] = [
    ("& (Suppli et al.,", "& Suppli et al.,"),
    ("& (Verma et al., 2024) &", "& Verma et al., 2024 &"),
    ("& (Kurki et al., 2023) &", "& Kurki et al., 2023 &"),
    ("& (Karczewski et al., 2025) &", "& Karczewski et al., 2025 &"),
    ("& (Boel et al., 2025) &", "& Boel et al., 2025 &"),
    ("& (Zhu et al., 2026) &", "& Zhu et al., 2026 &"),
    ("& (Dobin et al., 2013) &", "& Dobin et al., 2013 &"),
    ("& (Liao et al., 2014) &", "& Liao et al., 2014 &"),
    ("& (Robinson et al., 2010) &", "& Robinson et al., 2010 &"),
    ("& (Ritchie et al., 2015);", "& Ritchie et al., 2015;"),
    ("& (Stephens, 2017) &", "& Stephens, 2017 &"),
    ("& (Wang et al., 2020) &", "& Wang et al., 2020 &"),
    ("& (Dai et al., 2023) &", "& Dai et al., 2023 &"),
    ("& (Hinrichs et al., 2006; Lawrence et al., 2009) &", "& Hinrichs et al., 2006; Lawrence et al., 2009 &"),
    ("& (Oh et al., 2020) &", "& Oh et al., 2020 &"),
    ("& (DeTomaso and Yosef, 2021) &", "& DeTomaso and Yosef, 2021 &"),
    ("& (Yang et al., 2020) &", "& Yang et al., 2020 &"),
    ("& (Wolock et al., 2019) &", "& Wolock et al., 2019 &"),
    ("& (Badia-i-Mompel et al., 2022) &", "& Badia-i-Mompel et al., 2022 &"),
    ("& (Chu et al., 2022); (Hu and Chikina, 2024) &", "& Chu et al., 2022; Hu and Chikina, 2024 &"),
    ("& (Kleshchevnikov et al., 2022) &", "& Kleshchevnikov et al., 2022 &"),
    ("& (Traag et al., 2019); (McInnes et al., 2018) &", "& Traag et al., 2019; McInnes et al., 2018 &"),
    ("& (Hu et al., 2024); (Korsunsky et al., 2019) &", "& Hu et al., 2024; Korsunsky et al., 2019 &"),
    ("& (Domínguez Conde et al., 2022) &", "& Domínguez Conde et al., 2022 &"),
    ("& (Chang et al., 2015) &", "& Chang et al., 2015 &"),
    ("& (Zhang et al., 2024); (Zhang et al., 2008) &", "& Zhang et al., 2024; Zhang et al., 2008 &"),
    ("& (Sugimoto et al., 2025);", "& Sugimoto et al., 2025;"),
    ("(Broadaway et al., 2024) (Broadaway et al., 2024)", "(Broadaway et al., 2024)"),
]

# 4. Verified CrossRef/PubMed references to append to sections/references.tex (if DOI/substring not already present)
NEW_REFERENCES: list[tuple[str, str]] = [
    ("10.1152/ajpgi.00358.2018", "Suppli, M. P. et al.~Hepatic transcriptome signatures in patients with varying degrees of nonalcoholic fatty liver disease compared with healthy normal-weight individuals. \\emph{Am J Physiol Gastrointest Liver Physiol} \\textbf{316}, G462--G472 (2019). \\url{https://doi.org/10.1152/ajpgi.00358.2018}"),
    ("10.1016/j.jcmgh.2022.12.008", "Chen, T. et al.~Hepatocyte Smoothened activity controls susceptibility to insulin resistance and nonalcoholic fatty liver disease. \\emph{Cell Mol Gastroenterol Hepatol} \\textbf{15}, 949--970 (2023). \\url{https://doi.org/10.1016/j.jcmgh.2022.12.008}"),
    ("10.1126/scitranslmed.abo4474", "Fujiwara, N. et al.~Molecular signatures of long-term hepatocellular carcinoma risk in nonalcoholic fatty liver disease. \\emph{Sci Transl Med} \\textbf{14}, eabo4474 (2022). \\url{https://doi.org/10.1126/scitranslmed.abo4474}"),
    ("10.1097/HEP.0000000000001360", "Fujiwara, N. et al.~AI-based phenotyping of hepatic fiber morphology to inform molecular alterations in metabolic dysfunction-associated steatotic liver disease. \\emph{Hepatology} \\textbf{83}, 575--590 (2026). \\url{https://doi.org/10.1097/HEP.0000000000001360}"),
    ("10.7554/eLife.109534", "Zhao, Q. et al.~Global molecular landscape of early MASLD progression in human obesity. \\emph{eLife} \\textbf{14}, RP109534 (2026). \\url{https://doi.org/10.7554/eLife.109534}"),
    ("10.1038/s41467-023-37209-1", "Powell, N. R. et al.~Clinically important alterations in pharmacogene expression in histologically severe nonalcoholic fatty liver disease. \\emph{Nat Commun} \\textbf{14}, 1474 (2023). \\url{https://doi.org/10.1038/s41467-023-37209-1}"),
    ("10.1038/s41467-024-48956-0", "Verschuren, L. et al.~Development of a novel non-invasive biomarker panel for hepatic fibrosis in MASLD. \\emph{Nat Commun} \\textbf{15}, 4564 (2024). \\url{https://doi.org/10.1038/s41467-024-48956-0}"),
    ("10.1053/j.gastro.2015.05.039", "Xanthakos, S. A. et al.~High prevalence of nonalcoholic fatty liver disease in adolescents undergoing bariatric surgery. \\emph{Gastroenterology} \\textbf{149}, 623--634.e8 (2015). \\url{https://doi.org/10.1053/j.gastro.2015.05.039}"),
    ("10.1002/hep4.1940", "Yao, K. et al.~Transcriptomic profiling of a multiethnic pediatric NAFLD cohort reveals genes and pathways associated with disease. \\emph{Hepatol Commun} \\textbf{6}, 1598--1610 (2022). \\url{https://doi.org/10.1002/hep4.1940}"),
    ("10.1016/j.cmet.2013.07.004", "Ahrens, M. et al.~DNA methylation analysis in nonalcoholic fatty liver disease suggests distinct disease-specific and remodeling signatures after bariatric surgery. \\emph{Cell Metab} \\textbf{18}, 296--302 (2013). \\url{https://doi.org/10.1016/j.cmet.2013.07.004}"),
    ("10.1016/j.jhepr.2024.101167", "García-Nieto, E. et al.~Clinical and molecular characterization of steatotic liver disease in the setting of immune-mediated inflammatory diseases. \\emph{JHEP Rep} \\textbf{6}, 101167 (2024). \\url{https://doi.org/10.1016/j.jhepr.2024.101167}"),
    ("10.1016/j.xcrm.2025.101935", "Yang, M. et al.~Integrative proteo-transcriptomic characterization of advanced fibrosis in metabolic dysfunction-associated steatotic liver disease. \\emph{Cell Rep Med} \\textbf{6}, 101935 (2025). \\url{https://doi.org/10.1016/j.xcrm.2025.101935}"),
    ("10.1172/jci.insight.182232", "Steinhauser, M. L. et al.~The transcription factor ZNF469 regulates collagen production in liver fibrosis. \\emph{JCI Insight} \\textbf{10}, e182232 (2025). \\url{https://doi.org/10.1172/jci.insight.182232}"),
    ("10.1210/js.2018-00122", "Gerhard, G. S. et al.~Transcriptomic profiling of obesity-related nonalcoholic steatohepatitis (NASH). \\emph{J Endocr Soc} \\textbf{2}, 710--726 (2018). \\url{https://doi.org/10.1210/js.2018-00122}"),
    ("10.1126/science.aaz1776", "GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. \\emph{Science} \\textbf{369}, 1318--1330 (2020). \\url{https://doi.org/10.1126/science.aaz1776}"),
    ("10.1038/s41586-025-08677-w", "Sugimoto, A. et al.~Hepatic stellate cells control liver zonation, size and functions via R-spondin 3. \\emph{Nature} \\textbf{640}, 752--761 (2025). \\url{https://doi.org/10.1038/s41586-025-08677-w}"),
    ("10.1016/j.ajhg.2025.11.009", "Alkhawaja, A. A. et al.~Liver single-nucleus multiome profiling reveals cell-type mechanisms for cardiometabolic traits. \\emph{Am J Hum Genet} \\textbf{113}, 133--148 (2026). \\url{https://doi.org/10.1016/j.ajhg.2025.11.009}"),
    ("10.1038/s43856-025-00780-3", "Boel, F. et al.~Deep proteome profiling of metabolic dysfunction-associated steatotic liver disease. \\emph{Commun Med} \\textbf{5}, 56 (2025). \\url{https://doi.org/10.1038/s43856-025-00780-3}"),
    ("10.1016/j.cell.2025.11.031", "Tzouanas, C. N. et al.~Hepatic adaptation to chronic metabolic stress primes tumorigenesis. \\emph{Cell} \\textbf{189}, 435--460.e28 (2026). \\url{https://doi.org/10.1016/j.cell.2025.11.031}"),
    ("10.1126/science.adj1182", "Verma, A. et al.~Diversity and scale: genetic architecture of 2068 traits in the VA Million Veteran Program. \\emph{Science} \\textbf{385}, eadj1182 (2024). \\url{https://doi.org/10.1126/science.adj1182}"),
    ("10.1038/s41586-022-05473-8", "Kurki, M. I. et al.~FinnGen provides genetic insights from a well-phenotyped isolated population. \\emph{Nature} \\textbf{613}, 508--518 (2023). \\url{https://doi.org/10.1038/s41586-022-05473-8}"),
    ("10.1038/s41588-025-02335-7", "Karczewski, K. J. et al.~Pan-UK Biobank genome-wide association analyses enhance discovery and resolution of ancestry-enriched effects. \\emph{Nat Genet} \\textbf{57}, 2408--2417 (2025). \\url{https://doi.org/10.1038/s41588-025-02335-7}"),
    ("10.1038/s41588-024-01720-y", "Ghouse, J. et al.~Integrative common and rare variant analyses provide insights into the genetic architecture of liver cirrhosis. \\emph{Nat Genet} \\textbf{56}, 827--837 (2024). \\url{https://doi.org/10.1038/s41588-024-01720-y}"),
    ("10.1038/s41588-022-01270-1", "Chen, Y. et al.~Genomic atlas of the plasma metabolome prioritizes metabolites implicated in human diseases. \\emph{Nat Genet} \\textbf{55}, 44--53 (2023). \\url{https://doi.org/10.1038/s41588-022-01270-1}"),
    ("10.1038/s41467-023-42532-8", "Ottensmann, L. et al.~Genome-wide association analysis of plasma lipidome identifies 495 genetic associations. \\emph{Nat Commun} \\textbf{14}, 6934 (2023). \\url{https://doi.org/10.1038/s41467-023-42532-8}"),
    ("10.1038/s41590-023-01588-w", "Zhao, J. H. et al.~Genetics of circulating inflammatory proteins identifies drivers of immune-mediated disease risk and therapeutic targets. \\emph{Nat Immunol} \\textbf{24}, 1540--1551 (2023). \\url{https://doi.org/10.1038/s41590-023-01588-w}"),
    ("10.1038/s41467-025-57616-w", "Du, K. et al.~Targeting senescent hepatocytes for treatment of metabolic dysfunction-associated steatotic liver disease and multi-organ dysfunction. \\emph{Nat Commun} \\textbf{16}, 3038 (2025). \\url{https://doi.org/10.1038/s41467-025-57616-w}"),
    ("10.1172/JCI199303", "Qiu, X. et al.~THEMIS attenuates MASH by suppressing disease-associated hepatocyte induction and hepatocyte senescence in mice. \\emph{J Clin Invest} \\textbf{136}, e199303 (2026). \\url{https://doi.org/10.1172/JCI199303}"),
    ("10.1124/dmd.109.027466", "Fisher, C. D. et al.~Hepatic cytochrome P450 enzyme alterations in humans with progressive stages of nonalcoholic fatty liver disease. \\emph{Drug Metab Dispos} \\textbf{37}, 2087--2094 (2009). \\url{https://doi.org/10.1124/dmd.109.027466}"),
    ("10.1053/j.gastro.2013.07.047", "Murphy, S. K. et al.~Relationship between methylome and transcriptome in patients with nonalcoholic fatty liver disease. \\emph{Gastroenterology} \\textbf{145}, 1076--1087 (2013). \\url{https://doi.org/10.1053/j.gastro.2013.07.047}"),
    ("10.1016/j.jhep.2015.08.037", "Zubiete-Franco, I. et al.~Methionine and S-adenosylmethionine levels are critical regulators of PP2A activity modulating lipophagy during steatosis. \\emph{J Hepatol} \\textbf{64}, 409--418 (2016). \\url{https://doi.org/10.1016/j.jhep.2015.08.037}"),
    ("10.1016/j.molmet.2019.08.008", "Fernández-Tussy, P. et al.~miR-873-5p targets mitochondrial GNMT-Complex II interface contributing to non-alcoholic fatty liver disease. \\emph{Mol Metab} \\textbf{29}, 40--54 (2019). \\url{https://doi.org/10.1016/j.molmet.2019.08.008}"),
    ("10.1111/liv.14174", "Männistö, V. et al.~Total liver phosphatidylcholine content associates with non-alcoholic steatohepatitis and glycine N-methyltransferase expression. \\emph{Liver Int} \\textbf{39}, 1895--1905 (2019). \\url{https://doi.org/10.1111/liv.14174}"),
    ("10.1038/nbt.3269", "Gaidatzis, D., Burger, L., Florescu, M. \\& Stadler, M. B. Analysis of intronic and exonic reads in RNA-seq data characterizes transcriptional and post-transcriptional regulation. \\emph{Nat Biotechnol} \\textbf{33}, 722--729 (2015). \\url{https://doi.org/10.1038/nbt.3269}"),
    ("10.1038/nature15393", "1000 Genomes Project Consortium. A global reference for human genetic variation. \\emph{Nature} \\textbf{526}, 68--74 (2015). \\url{https://doi.org/10.1038/nature15393}"),
    ("10.1093/nar/gkac1071", "Frankish, A. et al.~GENCODE: reference annotation for the human and mouse genomes in 2023. \\emph{Nucleic Acids Res} \\textbf{51}, D942--D949 (2023). \\url{https://doi.org/10.1093/nar/gkac1071}"),
    ("10.1038/s41586-023-06045-0", "Chen, S. et al.~A genomic mutational constraint map using variation in 76,156 human genomes. \\emph{Nature} \\textbf{625}, 92--100 (2024). \\url{https://doi.org/10.1038/s41586-023-06045-0}"),
    ("10.1016/j.cell.2017.06.010", "Tsherniak, A. et al.~Defining a cancer dependency map. \\emph{Cell} \\textbf{170}, 564--576.e16 (2017). \\url{https://doi.org/10.1016/j.cell.2017.06.010}"),
    ("10.1126/science.1260419", "Uhlén, M. et al.~Tissue-based map of the human proteome. \\emph{Science} \\textbf{347}, 1260419 (2015). \\url{https://doi.org/10.1126/science.1260419}"),
    ("10.1534/g3.117.041277", "Hart, T. et al.~Evaluation and design of genome-wide CRISPR/SpCas9 knockout screens. \\emph{G3} \\textbf{7}, 2719--2727 (2017). \\url{https://doi.org/10.1534/g3.117.041277}"),
    ("10.1101/gr.240663.118", "Garcia-Alonso, L., Holland, C. H., Ibrahim, M. M., Turei, D. \\& Saez-Rodriguez, J. Benchmark and integration of resources for the estimation of human transcription factor activities. \\emph{Genome Res} \\textbf{29}, 1363--1375 (2019). \\url{https://doi.org/10.1101/gr.240663.118}"),
    ("10.1016/j.cels.2015.12.004", "Liberzon, A. et al.~The Molecular Signatures Database hallmark gene set collection. \\emph{Cell Syst} \\textbf{1}, 417--425 (2015). \\url{https://doi.org/10.1016/j.cels.2015.12.004}"),
    ("10.1093/bioinformatics/bts635", "Dobin, A. et al.~STAR: ultrafast universal RNA-seq aligner. \\emph{Bioinformatics} \\textbf{29}, 15--21 (2013). \\url{https://doi.org/10.1093/bioinformatics/bts635}"),
    ("10.1093/bioinformatics/btt656", "Liao, Y., Smyth, G. K. \\& Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. \\emph{Bioinformatics} \\textbf{30}, 923--930 (2014). \\url{https://doi.org/10.1093/bioinformatics/btt656}"),
    ("10.1093/bioinformatics/btp616", "Robinson, M. D., McCarthy, D. J. \\& Smyth, G. K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. \\emph{Bioinformatics} \\textbf{26}, 139--140 (2010). \\url{https://doi.org/10.1093/bioinformatics/btp616}"),
    ("10.1093/nar/gkv007", "Ritchie, M. E. et al.~limma powers differential expression analyses for RNA-sequencing and microarray studies. \\emph{Nucleic Acids Res} \\textbf{43}, e47 (2015). \\url{https://doi.org/10.1093/nar/gkv007}"),
    ("10.1186/gb-2014-15-2-r29", "Law, C. W., Chen, Y., Shi, W. \\& Smyth, G. K. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. \\emph{Genome Biol} \\textbf{15}, R29 (2014). \\url{https://doi.org/10.1186/gb-2014-15-2-r29}"),
    ("10.1093/nar/gkv412", "Liu, R. et al.~Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. \\emph{Nucleic Acids Res} \\textbf{43}, e97 (2015). \\url{https://doi.org/10.1093/nar/gkv412}"),
    ("10.1093/biostatistics/kxw041", "Stephens, M. False discovery rates: a new deal. \\emph{Biostatistics} \\textbf{18}, 275--294 (2017). \\url{https://doi.org/10.1093/biostatistics/kxw041}"),
    ("10.1111/rssb.12388", "Wang, G., Sarkar, A., Carbonetto, P. \\& Stephens, M. A simple new approach to variable selection in regression, with application to genetic fine mapping. \\emph{J R Stat Soc Series B Stat Methodol} \\textbf{82}, 1273--1300 (2020). \\url{https://doi.org/10.1111/rssb.12388}"),
    ("10.1371/journal.pgen.1004383", "Giambartolomei, C. et al.~Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. \\emph{PLoS Genet} \\textbf{10}, e1004383 (2014). \\url{https://doi.org/10.1371/journal.pgen.1004383}"),
    ("10.1371/journal.pgen.1009440", "Wallace, C. A more accurate method for colocalisation analysis allowing for multiple causal variants. \\emph{PLoS Genet} \\textbf{17}, e1009440 (2021). \\url{https://doi.org/10.1371/journal.pgen.1009440}"),
    ("10.1038/s41588-020-00735-5", "Weissbrod, O. et al.~Functionally informed fine-mapping and polygenic localization of complex trait heritability. \\emph{Nat Genet} \\textbf{52}, 1355--1363 (2020). \\url{https://doi.org/10.1038/s41588-020-00735-5}"),
    ("10.1038/s41467-023-36862-w", "Dai, Q. et al.~OTTERS: a powerful TWAS framework leveraging summary-level reference data. \\emph{Nat Commun} \\textbf{14}, 1271 (2023). \\url{https://doi.org/10.1038/s41467-023-36862-w}"),
    ("10.1038/s41588-024-01870-z", "Yuan, K. et al.~Fine-mapping across diverse ancestries drives the discovery of putative causal variants underlying human complex traits and diseases. \\emph{Nat Genet} \\textbf{56}, 1841--1850 (2024). \\url{https://doi.org/10.1038/s41588-024-01870-z}"),
    ("10.1038/s41588-023-01604-7", "Gao, B. \\& Zhou, X. MESuSiE enables scalable and powerful multi-ancestry fine-mapping of causal variants in genome-wide association studies. \\emph{Nat Genet} \\textbf{56}, 170--179 (2024). \\url{https://doi.org/10.1038/s41588-023-01604-7}"),
    ("10.1186/s13059-016-0974-4", "McLaren, W. et al.~The Ensembl Variant Effect Predictor. \\emph{Genome Biol} \\textbf{17}, 122 (2016). \\url{https://doi.org/10.1186/s13059-016-0974-4}"),
    ("10.1093/nar/gkj144", "Hinrichs, A. S. et al.~The UCSC Genome Browser Database: update 2006. \\emph{Nucleic Acids Res} \\textbf{34}, D590--D598 (2006). \\url{https://doi.org/10.1093/nar/gkj144}"),
    ("10.1093/bioinformatics/btp328", "Lawrence, M., Gentleman, R. \\& Carey, V. rtracklayer: an R package for interfacing with genome browsers. \\emph{Bioinformatics} \\textbf{25}, 1841--1842 (2009). \\url{https://doi.org/10.1093/bioinformatics/btp328}"),
    ("10.12688/f1000research.28033.1", "Oh, S. et al.~HGNChelper: identification and correction of invalid gene symbols for human and mouse. \\emph{F1000Res} \\textbf{9}, 1493 (2020). \\url{https://doi.org/10.12688/f1000research.28033.1}"),
    ("10.1016/j.cels.2021.04.005", "DeTomaso, D. \\& Yosef, N. Hotspot identifies informative gene modules across modalities of single-cell genomics. \\emph{Cell Syst} \\textbf{12}, 446--456.e9 (2021). \\url{https://doi.org/10.1016/j.cels.2021.04.005}"),
    ("10.1038/s41592-018-0229-2", "Lopez, R., Regier, J., Cole, M. B., Jordan, M. I. \\& Yosef, N. Deep generative modeling for single-cell transcriptomics. \\emph{Nat Methods} \\textbf{15}, 1053--1058 (2018). \\url{https://doi.org/10.1038/s41592-018-0229-2}"),
    ("10.1186/s13059-020-1950-6", "Yang, S. et al.~Decontamination of ambient RNA in single-cell RNA-seq with DecontX. \\emph{Genome Biol} \\textbf{21}, 57 (2020). \\url{https://doi.org/10.1186/s13059-020-1950-6}"),
    ("10.1016/j.cels.2018.11.005", "Wolock, S. L., Lopez, R. \\& Klein, A. M. Scrublet: computational identification of cell doublets in single-cell transcriptomic data. \\emph{Cell Syst} \\textbf{8}, 281--291.e9 (2019). \\url{https://doi.org/10.1016/j.cels.2018.11.005}"),
    ("10.1038/s41467-022-30755-0", "Dimitrov, D. et al.~Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. \\emph{Nat Commun} \\textbf{13}, 3224 (2022). \\url{https://doi.org/10.1038/s41467-022-30755-0}"),
    ("10.1093/bioadv/vbac016", "Badia-i-Mompel, P. et al.~decoupleR: ensemble of computational methods to infer biological activities from omics data. \\emph{Bioinform Adv} \\textbf{2}, vbac016 (2022). \\url{https://doi.org/10.1093/bioadv/vbac016}"),
    ("10.1038/s43018-022-00356-3", "Chu, T., Wang, Z., Pe'er, D. \\& Danko, C. G. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. \\emph{Nat Cancer} \\textbf{3}, 505--517 (2022). \\url{https://doi.org/10.1038/s43018-022-00356-3}"),
    ("10.1093/bioinformatics/btae440", "Hu, M. \\& Chikina, M. InstaPrism: an R package for fast implementation of BayesPrism. \\emph{Bioinformatics} \\textbf{40}, btae440 (2024). \\url{https://doi.org/10.1093/bioinformatics/btae440}"),
    ("10.1038/s41587-021-01139-4", "Kleshchevnikov, V. et al.~Cell2location maps fine-grained cell types in spatial transcriptomics. \\emph{Nat Biotechnol} \\textbf{40}, 661--671 (2022). \\url{https://doi.org/10.1038/s41587-021-01139-4}"),
    ("10.1038/s41598-019-41695-z", "Traag, V. A., Waltman, L. \\& van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. \\emph{Sci Rep} \\textbf{9}, 5233 (2019). \\url{https://doi.org/10.1038/s41598-019-41695-z}"),
    ("10.21105/joss.00861", "McInnes, L., Healy, J., Saul, N. \\& Großberger, L. UMAP: uniform manifold approximation and projection. \\emph{J Open Source Softw} \\textbf{3}, 861 (2018). \\url{https://doi.org/10.21105/joss.00861}"),
    ("10.1093/bioadv/vbaf167", "Hu, M. et al.~ScaleSC: a superfast and scalable single-cell RNA-seq analysis tool. \\emph{Bioinform Adv} \\textbf{5}, vbaf167 (2024). \\url{https://doi.org/10.1093/bioadv/vbaf167}"),
    ("10.1038/s41592-019-0619-0", "Korsunsky, I. et al.~Fast, sensitive and accurate integration of single-cell data with Harmony. \\emph{Nat Methods} \\textbf{16}, 1289--1296 (2019). \\url{https://doi.org/10.1038/s41592-019-0619-0}"),
    ("10.1126/science.abl5197", "Domínguez Conde, C. et al.~Cross-tissue immune cell analysis reveals tissue-specific features in humans. \\emph{Science} \\textbf{376}, eabl5197 (2022). \\url{https://doi.org/10.1126/science.abl5197}"),
    ("10.1186/s13742-015-0047-8", "Chang, C. C. et al.~Second-generation PLINK: rising to the challenge of larger and richer datasets. \\emph{GigaScience} \\textbf{4}, s13742-015-0047-8 (2015). \\url{https://doi.org/10.1186/s13742-015-0047-8}"),
    ("10.1038/s41592-023-02139-9", "Zhang, K. et al.~A fast, scalable and versatile tool for analysis of single-cell omics data. \\emph{Nat Methods} \\textbf{21}, 217--227 (2024). \\url{https://doi.org/10.1038/s41592-023-02139-9}"),
    ("10.1186/gb-2008-9-9-r137", "Zhang, Y. et al.~Model-based analysis of ChIP-Seq (MACS). \\emph{Genome Biol} \\textbf{9}, R137 (2008). \\url{https://doi.org/10.1186/gb-2008-9-9-r137}"),
    ("10.1038/s41562-019-0566-x", "Grotzinger, A. D. et al.~Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. \\emph{Nat Hum Behav} \\textbf{3}, 513--525 (2019). \\url{https://doi.org/10.1038/s41562-019-0566-x}"),
    ("10.1126/science.aec8514", "Pearce, J. D. et al.~TranscriptFormer: a generative cell atlas across 1.5 billion years of evolution. \\emph{Science} \\textbf{393}, aec8514 (2026). \\url{https://doi.org/10.1126/science.aec8514}"),
    ("10.1101/2024.08.16.608180", "Chen, H. \\& Theodoris, C. V. Quantized multi-task learning for context-specific representations of gene network dynamics. Preprint at bioRxiv \\url{https://doi.org/10.1101/2024.08.16.608180} (2024)."),
    ("10.1101/2025.06.11.659222", "Kang, B. et al.~A large-scale foundation model for bulk transcriptomes. Preprint at bioRxiv \\url{https://doi.org/10.1101/2025.06.11.659222} (2025)."),
    ("arXiv:2306.15794", "Nguyen, E. et al.~HyenaDNA: long-range genomic sequence modeling at single nucleotide resolution. \\emph{Adv Neural Inf Process Syst} \\textbf{36}, 43177--43201 (2023). \\url{https://arxiv.org/abs/2306.15794}"),
    ("10.1038/s41592-024-02523-z", "Dalla-Torre, H. et al.~Nucleotide Transformer: building and evaluating robust foundation models for human genomics. \\emph{Nat Methods} \\textbf{22}, 287--297 (2025). \\url{https://doi.org/10.1038/s41592-024-02523-z}"),
    ("10.1101/2023.11.28.568918", "Rosen, Y. et al.~Universal Cell Embeddings: a foundation model for cell biology. Preprint at bioRxiv \\url{https://doi.org/10.1101/2023.11.28.568918} (2023)."),
    ("10.1038/s41586-024-08411-y", "Heimberg, G. et al.~A cell atlas foundation model for scalable search of similar human cells. \\emph{Nature} \\textbf{638}, 1085--1094 (2025). \\url{https://doi.org/10.1038/s41586-024-08411-y}"),
    ("10.1080/07350015.1983.10509354", "Freedman, D. \\& Lane, D. A nonstochastic interpretation of reported significance levels. \\emph{J Bus Econ Stat} \\textbf{1}, 292--298 (1983). \\url{https://doi.org/10.1080/07350015.1983.10509354}"),
    ("10.1080/00031305.2000.10474549", "Long, J. S. \\& Ervin, L. H. Using heteroscedasticity consistent standard errors in the linear regression model. \\emph{Am Stat} \\textbf{54}, 217--224 (2000). \\url{https://doi.org/10.1080/00031305.2000.10474549}"),
    ("10.1002/sim.1482", "Knapp, G. \\& Hartung, J. Improved tests for a random effects meta-regression with a single covariate. \\emph{Stat Med} \\textbf{22}, 2693--2710 (2003). \\url{https://doi.org/10.1002/sim.1482}"),
    ("10.1093/jnci/22.4.719", "Mantel, N. \\& Haenszel, W. Statistical aspects of the analysis of data from retrospective studies of disease. \\emph{J Natl Cancer Inst} \\textbf{22}, 719--748 (1959). \\url{https://doi.org/10.1093/jnci/22.4.719}"),
    ("10.2307/2531052", "Robins, J., Breslow, N. \\& Greenland, S. Estimators of the Mantel-Haenszel variance consistent in both sparse data and large-strata limiting models. \\emph{Biometrics} \\textbf{42}, 311--323 (1986). \\url{https://doi.org/10.2307/2531052}"),
    ("10.1093/biomet/37.1-2.17", "Moran, P. A. P. Notes on continuous stochastic phenomena. \\emph{Biometrika} \\textbf{37}, 17--23 (1950). \\url{https://doi.org/10.1093/biomet/37.1-2.17}"),
]

PLACEHOLDER_RE = re.compile(r"\{\[\}[^}]*?\{\]\}")
NUMERIC_INTERVAL_RE = re.compile(r"^\{\[\}[-+\u22120-9.,;\s–—~%]+?\{\]\}$")


def find_unresolved(text: str) -> list[str]:
    out = []
    for m in PLACEHOLDER_RE.finditer(text):
        tok = m.group(0)
        if not NUMERIC_INTERVAL_RE.match(tok):
            out.append(tok)
    return out


def transform_text(rel_path: str, text: str) -> tuple[str, int]:
    new_text = text
    n_ops = 0
    for old, new in GLOBAL_REPLACEMENTS:
        if old in new_text:
            c = new_text.count(old)
            new_text = new_text.replace(old, new)
            n_ops += c
    for old, new in FILE_SPECIFIC_REPLACEMENTS.get(rel_path, []):
        if old in new_text:
            c = new_text.count(old)
            new_text = new_text.replace(old, new)
            n_ops += c
    for old, new in KRT_CLEANUPS:
        if old in new_text:
            new_text = new_text.replace(old, new)
    return new_text, n_ops


def main() -> int:
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    tex_files = sorted(TEX_DIR.rglob("*.tex"))
    total_before = 0
    total_after = 0
    unresolved_all: list[tuple[str, str]] = []
    changed_files: list[tuple[str, int]] = []

    for p in tex_files:
        rel = str(p.relative_to(ROOT))
        orig = p.read_text(encoding="utf-8")
        before = find_unresolved(orig)
        total_before += len(before)
        updated, n_ops = transform_text(rel, orig)
        after = find_unresolved(updated)
        total_after += len(after)
        for tok in after:
            unresolved_all.append((rel, tok))
        if updated != orig:
            changed_files.append((rel, n_ops))
            if args.apply:
                p.write_text(updated, encoding="utf-8")

    # References update
    refs_text = REFS_FILE.read_text(encoding="utf-8")
    to_add = [entry for key, entry in NEW_REFERENCES if key not in refs_text]
    if to_add and args.apply:
        assert "\\end{enumerate}" in refs_text
        block = "\n".join(f"\\item\n  {e}" for e in to_add) + "\n\\end{enumerate}"
        REFS_FILE.write_text(refs_text.replace("\\end{enumerate}", block, 1), encoding="utf-8")

    print(f"Mode: {'APPLY' if args.apply else 'DRY-RUN'}")
    print(f"Scanned {len(tex_files)} .tex files")
    print(f"Non-interval placeholders: {total_before} -> {total_after}")
    print(f"Files modified: {len(changed_files)}")
    for rel, n in changed_files:
        print(f"  {rel}: {n} replacements")
    print(f"New verified references to append to references.tex: {len(to_add)}")
    if unresolved_all:
        print("UNRESOLVED PLACEHOLDERS:")
        for rel, tok in unresolved_all:
            print(f"  {rel}: {tok}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
