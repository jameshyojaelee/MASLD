"""Build the 'cat' text-writer edit list (Gene Catalog text, framing, references).

Checks every old string occurs exactly once in its source file and that no two
edits overlap in the same file. Writes RUN/staged/edits_cat.json and, after
simulating the edits in memory, RUN/staged/ref_placeholders.tsv.
Read-only on docs/.
"""
import json, re, sys
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z"
R9D = ("Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
       "catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet")
C12 = "RUN/inputs/C12-catalog.json"
C14 = "RUN/inputs/C14-framing.json"
M = "docs/manuscript/METHODS.md"
L = "docs/manuscript/05_figure_legends.md"
D = "docs/manuscript/draft/"

E = []
def add(id_, file, old, new, items, reason, number_sources=None, decision=None):
    E.append(dict(id=id_, file=file, old=old, new=new, items=items,
                  number_sources=number_sources or [], decision=decision, reason=reason))

# ---------------- A6.1: "unsupported" = tested without support ----------------
add("cat-01", D + "abstract.md",
    "distinguishes untested from unsupported, and assigns each gene a next experiment.",
    "separates untested genes from genes tested without support, and proposes a next experiment for each gene; for 94.5% of genes, that experiment is to obtain the missing evidence layer.",
    ["A6.1", "A6.2"],
    "'Unsupported' can be read as a negative call; most routes only name the missing layer.",
    [f"{R9D}: next_experiment_rule_id EXP_DIRECT_ASSAY_V1 13,551 + EXP_UNRESOLVED_V1 13,127 = 26,678 of 28,233 (94.49%)"])
add("cat-02", D + "abstract.md",
    "Most disease-state calls without colocalization were untested, not genetic negatives.",
    "Most disease-state calls without colocalization had no multi-signal colocalization result, so they are not genetic negatives.",
    ["B-C2"],
    "The Catalog gives these 935 genes an indeterminate (tested, ABF-only) genetic call, not untested.")
add("cat-03", D + "01_intro.md",
    "and assigns the next discriminating experiment (Figure 6).",
    "and proposes a next experiment (Figure 6). For 26,678 of the 28,233 genes (94.5%), that experiment is to obtain the missing evidence layer.",
    ["A6.2"],
    "States how many routes are generic.",
    [f"{R9D}: next_experiment_rule_id EXP_DIRECT_ASSAY_V1 13,551 + EXP_UNRESOLVED_V1 13,127"])
add("cat-04", D + "01_intro.md",
    "separates untested from unsupported, is the standalone contribution",
    "separates untested genes from genes tested without support, is the standalone contribution",
    ["A6.1"], "Same wording as the abstract.")
add("cat-05", D + "fig5_discussion.md",
    "No nonsignificant result is converted into a negative call.",
    "No nonsignificant result is converted into a negative call. In this paper,\nunsupported means tested without support (indeterminate), not tested negative.",
    ["A6.1"], "Defines 'unsupported' once; no negative-call rule is introduced (D7).")

# ---------------- B-C2: r8 label, 15,009 / 290 / 18, GLP1R ----------------
add("cat-06", D + "fig5_discussion.md",
    "The 9,312 genes with no Tier-1/2 COLOC record carry an underpowered-source label.\nNo power assessment was made for them, so the label marks a missing record, not a\nmeasured lack of power.",
    "The 9,312 genes with no Tier-1/2 COLOC record are genetics-untestable, with the\nreason recorded as no Tier-1/2 COLOC record. No power assessment was made for\nthem, so this state marks a missing record, not a measured lack of power.",
    ["B-C2"], "Removes the r8 'underpowered-source' label; matches r9d gen_untestable_reason_r9 and the Fig 6B legend.")
add("cat-07", D + "fig5_discussion.md",
    "best PP.H4 at most 0.5, the pattern expected when the trait and expression\nsignals have distinct causal variants. Eighteen of these 20 are called",
    "best PP.H4 at most 0.5. All 20 had SuSiE PP.H4 above 0.5 in at least one\nTier-1/2 study before the check removed signal pairs with too little shared\nposterior, so their PP.H3 pattern follows from which pairs the check kept and is\nnot independent evidence for distinct causal variants. Eighteen of these 20 are called",
    ["B-C2", "B-G1"],
    "The 20 'distinct-variant' genes were SuSiE-positive until the eligibility check (PLAN §1 row 3).",
    ["GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/qc/call_changes.tsv: 20 of 20 genes r2_status eligible_pp4_le_0.5, minimum earlier PP.H4 0.511 (C12 new_issues[1]; critic B-C2)"])
add("cat-08", D + "fig5_discussion.md",
    "Of these 3,335,\n3,319 are called indeterminate.",
    "Of these 3,335,\n3,319 are called indeterminate. The indeterminate genetic call also covers\n15,009 genes with a Tier-1/2 COLOC record but no SuSiE result; for them the\nmulti-signal test is missing, not negative. Among these, 290 unresolved and 18\nestablished-state genes have a single-signal (ABF) PP.H4 above 0.5, TM6SF2\n(0.64) among them.",
    ["B-C2"],
    "Discloses that 'indeterminate' includes 15,009 genes never tested by SuSiE, 308 of them ABF-supported.",
    [f"{R9D}: gen_susie_state_t12 no_susie_result 15,009",
     f"{R9D}: primary_evidence_class x gen_susie_state_t12 where gen_coloc_best_abf_pp4 > 0.5: unresolved/no_susie_result 290, established_state_associated/no_susie_result 18",
     f"{R9D}: TM6SF2 gen_coloc_best_abf_pp4 0.639"])
add("cat-09", D + "fig5_discussion.md",
    "GLP1R is untested. Its RNA was not detected and it has no Tier-1/2\nSuSiE result, so it routes first to the missing decisive assay.",
    "GLP1R is untested because its RNA was not detected. Its genetic call is\nindeterminate, with no Tier-1/2 SuSiE result, and it routes first to the missing\ndecisive assay.",
    ["B-C2"], "GLP1R's genetic call is indeterminate (r9d no_susie_result), the same state as MTARC1's.")
add("cat-10", L,
    "*GLP1R* is untested: RNA not detected and no Tier-1/2 SuSiE result.",
    "*GLP1R* is untested because its RNA was not detected; its genetic call is indeterminate (no Tier-1/2 SuSiE result).",
    ["B-C2"], "Same fix in the Fig 6C legend.")
add("cat-11", L,
    "The 20 unresolved genes with a coding variant at PIP ≥ 0.5 are routed by their eQTL-colocalization class.",
    "The 20 unresolved genes with a coding variant at PIP ≥ 0.5 are routed by their eQTL-colocalization class. Of the unresolved and established-state genes, 290 and 18 have ABF PP.H4 > 0.5 but no Tier-1/2 SuSiE result.",
    ["B-C2"], "Per-class count of ABF-supported genes without a SuSiE test (critic C12 missed[0]).",
    [f"{R9D}: primary_evidence_class x gen_susie_state_t12 where gen_coloc_best_abf_pp4 > 0.5"])
add("cat-12", D + "supplementary_notes.md",
    "These genes are untested by multi-signal colocalization, not genetic negatives.",
    "These genes have no multi-signal colocalization result. Their Catalog genetic call is indeterminate, and they are not genetic negatives.",
    ["B-C2"], "Aligns with the Catalog label for the 935 genes.")
add("cat-13", "docs/PAPER.md",
    "have no\n   multi-signal COLOC result and are\n   untested, not genetic negatives.",
    "have no\n   multi-signal COLOC result, so their genetic call is indeterminate, not\n   negative.",
    ["B-C2"], "Same alignment in the paper authority.")

# ---------------- A6.2: routes ----------------
add("cat-14", D + "fig5_discussion.md",
    "model. Unobservable candidates first need the missing decisive assay. Transcript",
    "model. Untested candidates (13,551) route to the missing decisive assay and\nunresolved candidates (13,127) to measurement of the most decisive missing layer;\nfor these 26,678 genes the route names the missing layer, not a specific assay. Transcript",
    ["A6.2"], "Says how many routes are generic; the unresolved route was not described.",
    [f"{R9D}: next_experiment_rule_id EXP_DIRECT_ASSAY_V1 13,551, EXP_UNRESOLVED_V1 13,127"])
add("cat-15", D + "fig5_discussion.md",
    "context-specific challenge remain separate routes.",
    "context-specific challenge are separate designs that no current record is\nassigned.",
    ["A6.2"], "No record carries these routes (molecular_object_sidecar.jsonl).")
add("cat-16", D + "fig5_discussion.md",
    "No gene is\ntherefore assigned to allele-aware testing of a regulatory element, although\nthat route is defined.",
    "The 400 genetically\nanchored genes route to allele-aware regulatory editing at the gene level, but no\nrecord names the element to edit, and the element-level route defined in Methods\nis assigned to no record. For 108 of the 426 genes with genetic support, the caQTL\nfield gives an accessible region that also colocalizes with the gene's eQTL.",
    ["A6.2"], "Removes the contradiction with lines 117-118.",
    [f"{R9D}: next_experiment_rule_id EXP_ALLELE_AWARE_V1 400; caqtl_gene_link M1_three_way_supported_linked 108 = 101 anchored + 7 concordant (C12 A6.2)",
     "426 = 400 genetically anchored + 26 concordant (fig5_discussion.md:116-117)"])
add("cat-17", D + "fig5_discussion.md",
    "Untestable or source-dependent candidates first need the missing decisive\nmeasurement.",
    "Untested and unresolved candidates, 26,678 of the 28,233 genes, first need the\nmissing decisive measurement.",
    ["A6.2"], "The closing paragraph implied tailored routes for all genes.",
    [f"{R9D}: next_experiment_rule_id EXP_DIRECT_ASSAY_V1 + EXP_UNRESOLVED_V1"])
add("cat-18", M,
    "Protein-coding genes keep the activation, inhibition, protein-level, lineage-restricted and context-specific routes.",
    "No current record is a regulatory element, so the regulatory-element route is assigned to none. Every protein-coding molecular-object record (14,176) receives one route, gene perturbation with rescue; activation, protein-level, lineage-restricted and context-specific designs are not assigned separately.",
    ["A6.2"], "All 14,176 protein-coding objects carry EXP_PROTEIN_STATE_CONTEXT_V2; no DNA-element record exists.",
    ["figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/molecular_object_sidecar.jsonl: candidate_object protein_coding_gene 14,176 (C12 A6.2; critic)"])
add("cat-19", L,
    "plus one molecular object (gene, RNA or DNA element) and one proposed discriminating experiment.",
    "plus one molecular object (protein-coding gene, lncRNA transcript, ambiguous lncRNA locus or untestable nomination; no record is a DNA element) and one proposed discriminating experiment.",
    ["A6.2"], "No record is a DNA element.")

# ---------------- A6.3: priors and p12 sensitivity; RNA threshold sensitivity ----------------
add("cat-20", M,
    "The prespecified high-resolution support threshold was SuSiE PP.H4 > 0.5.",
    "Both methods used the coloc priors p1 = p2 = 1 × 10⁻⁴ and p12 = 5 × 10⁻⁶. The prespecified high-resolution support threshold was SuSiE PP.H4 > 0.5. Because p12 enters only the shared-variant hypothesis, PP.H4 at another p12 follows exactly from the stored posteriors as k·PP.H4/(k·PP.H4 + 1 − PP.H4), where k is the ratio of the new to the original p12. At p12 = 1 × 10⁻⁶, 122 of the 436 supported genes fall to PP.H4 ≤ 0.5, and at p12 = 1 × 10⁻⁵, 29 genes tested below 0.5 rise above it (Supplementary Note 1). No Catalog class was changed.",
    ["A6.3", "B-C4"], "Priors were stated only in the caQTL paragraph; the p12 sensitivity for the 436 was unreported.",
    ["GWAS/finemapping/src/06_susie_coloc.R:111-113 (1e-4, 1e-4, 5e-6); coloc 5.2.3 defaults (C12 critic A6.3)",
     f"{R9D}: gen_coloc_best_susie_pp4 reweighted, k = 0.2: 122 of 436 fall to <= 0.5; k = 2: 29 tested <= 0.5 rise above 0.5 (recomputed by cat)"])
add("cat-21", D + "supplementary_notes.md",
    "Of the 436 primary multi-signal colocalized genes, 423 take their strongest colocalization evidence from a liver-enzyme trait and 13 from a direct disease or liver-fat trait.",
    "Of the 436 primary multi-signal colocalized genes, 423 take their strongest colocalization evidence from a liver-enzyme trait and 13 from a direct disease or liver-fat trait. Membership of genes near the threshold depends on the colocalization prior p12. At p12 = 1 × 10⁻⁶ instead of the prespecified 5 × 10⁻⁶, 122 of the 436 genes, those with best PP.H4 between 0.50 and 0.83, fall to PP.H4 ≤ 0.5: 111 of the 400 genetically anchored genes, 8 of the 26 supported on both axes and 3 untested genes with no RNA call. At p12 = 1 × 10⁻⁵, 29 unresolved genes rise above 0.5. The Catalog classes use the prespecified prior, and no class was changed.",
    ["A6.3"], "Reports the p12 sensitivity for gene membership; 5e-6 stays primary.",
    [f"{R9D}: k = 0.2 losers by primary_evidence_class genetically_anchored 111, concordant 8, untested 3; best PP.H4 range 0.5037-0.8310; k = 2 gainers unresolved 29 (recomputed by cat)"])
add("cat-22", D + "supplementary_notes.md",
    "losses concentrated in calls near the 0.5 log2FC threshold.",
    "losses concentrated in calls near the 0.5 log2FC threshold. The call count depends more on the fold-change floor than on the q threshold. In the Gene Catalog gene table, which holds 1,853 of the calls, a floor of 0.4 or 0.6 log2FC at BH q < 0.05 gives 2,826 or 1,229 calls, whereas 22 calls have q between 0.01 and 0.05.",
    ["A6.3"], "RNA threshold sensitivity, computed on the same fit (no refit).",
    [f"{R9D}: rna_bh_q < 0.05 & |rna_log2fc| > 0.5: 1,853; > 0.4: 2,826; > 0.6: 1,229; 0.01 <= q < 0.05 & |log2FC| > 0.5: 22 (recomputed by cat)"])

# ---------------- B-C4: SE in gene table ----------------
add("cat-23", M,
    "The evidence rows carry the native estimate and, for RNA, the BH q, but no standard error, interval, p value or unit count.",
    "The evidence rows carry the native estimate and, for RNA, the BH q, but no standard error, interval, p value or unit count. The gene table adds, for RNA, the limma standard error, moderated t statistic and p value, and the unit counts of the fit (774 genetic-identity clusters: 138 control, 636 disease). PP.H4 is a posterior probability and has no standard error.",
    ["B-C4"], "The r9d gene table already holds rna_se, rna_t, rna_p and fit-level unit counts.",
    [f"{R9D}: rna_n_units 774, rna_n_control_units 138, rna_n_disease_units 636 (C12 B-C4)"])

# ---------------- B-F1: KRT identifiers ----------------
add("cat-24", M,
    "| Repeat-biopsy pairs (source-labelled; identity not verified) | [ref: GSE193066 source study] | GEO: GSE193066 |",
    "| Disease-only cohort counted in the nine-cohort library census | Kozumi et al., 2021 | GEO: GSE167523 |\n| Repeat-biopsy pairs (source-labelled; identity not verified) | [ref: GSE193066 source study] | GEO: GSE193066 |",
    ["B-F1"], "GSE167523 is one of the nine census cohorts (fig1.md:5) but had no KRT row. GEO lists PMID 34105780 (Kozumi 2021); 98 biopsy-proven NAFLD patients, no controls.",
    ["GEO GSE167523 record (wget 2026-10-21): Series_pubmed_id 34105780"])
add("cat-25", M,
    "| **Deposited data: single-cell, chromatin and spatial** | | |",
    "| **Deposited data: liver expression microarray** | | |\n| Treatment-response cohort, paired biopsies (Figure 7E) | Lefebvre et al., 2017 | GEO: GSE83452 (Affymetrix Human Gene 2.0 ST Array, GPL16686) |\n| **Deposited data: single-cell, chromatin and spatial** | | |",
    ["B-F1"], "GSE83452 supplies Fig 7E but had no KRT row; it is a microarray, so it is not filed under bulk RNA-seq (critic C14).",
    ["GEO GPL16686 record (wget 2026-10-21): '[HuGene-2_0-st] Affymetrix Human Gene 2.0 ST Array'", "GEO GSE83452 record: Series_pubmed_id 28679947"])
add("cat-26", M,
    "Gribben et al., 2024 (GSE202379) [ref]; Andrews et al., 2022 (GSE185477) [ref]; Filliol et al. (GSE174748) [ref]; Guilliams et al., 2022 (GSE192742) [ref];",
    "Gribben et al., 2024 (GSE202379); Andrews et al., 2022 (GSE185477); Filliol et al., 2022 (GSE174748); Guilliams et al., 2022 (GSE192742); Kim et al., 2025 (GSE244832);",
    ["B-F3", "B4"], "Refs 28, 29, 15, 16 exist; adds the GSE244832 source study (GEO PMID 39522884).",
    ["GEO GSE244832 record: Series_pubmed_id 39522884"])
add("cat-27", M,
    "| Independent single-nucleus RNA-seq | [ref: source studies] | GEO: GSE256398, GSE289173 |",
    "| Independent single-nucleus RNA-seq | [ref: GSE256398 source study]; Hong et al., 2025 (GSE289173) | GEO: GSE256398, GSE289173 |",
    ["B4"], "GSE289173 source study (GEO PMID 40562914).", ["GEO GSE289173 record: Series_pubmed_id 40562914"])
add("cat-28", M,
    "| Single-nucleus ATAC-seq | [ref: source studies] | GEO: GSE244832, GSE281367 |",
    "| Single-nucleus ATAC-seq | Kim et al., 2025 (GSE244832); [ref: GSE281367 source study] | GEO: GSE244832, GSE281367 |",
    ["B4"], "GSE244832 source study.")
add("cat-29", M,
    "| Visium and spatial transcriptomics | Guilliams et al. (GSE192741) [ref]; Vu et al. [ref]; HMSMA [ref]; Yakubovsky et al. (2026) [ref] | GEO: GSE192741, GSE287826; GSA-Human: HRA007511 |",
    "| Visium and spatial transcriptomics | Guilliams et al., 2022 (GSE192741); Vu et al., 2026; Li et al., 2025 (HMSMA); Yakubovsky et al., 2026 | GEO: GSE192741; UQ eSpace: doi:10.48610/e95155f (Vu et al., CytAssist FFPE Visium); Zenodo: doi:10.5281/zenodo.17735506 and doi:10.5281/zenodo.17735587 (Yakubovsky et al.); GSA-Human: HRA007511 (HMSMA) |\n| GeoMx hepatocyte profiles (obtained, not analyzed) | [ref: GSE287826 source study] | GEO: GSE287826 |",
    ["B-F1"], "Identifiers for the analyzed Visium sources; GSE287826 marked as not analyzed (METHODS.md:419, :427).",
    ["Analysis/Spatial/config/spatial_datasets.yaml:152 (doi.org/10.48610/e95155f; resolves to espace.library.uq.edu.au/view/UQ:e95155f, wget 2026-10-21)",
     "Zenodo API records 17735506 and 17735587 (wget 2026-10-21): Yakubovsky and Itzkovitz, live-donor liver spatial atlas"])
add("cat-30", M,
    "| CosMx panel | Boesch et al., 2026 [ref] |",
    "| CosMx panel | Boesch et al., 2026 |", ["B-F3"], "Ref 18; duplicate marker.")
add("cat-31", M,
    "| Published liver-trait GWAS | Source studies |",
    "| Published liver-trait GWAS | Namjou et al., 2019; Anstee et al., 2020; Liu et al., 2021; Ghodsian et al., 2021; Haas et al., 2021; van der Meer et al., 2022; Sveinbjornsson et al., 2022 |",
    ["B-F3"], "Cites refs 30, 31, 32, 11, 33, 34 (otherwise uncited) and 26 by name; PMIDs map in the same order (critic C14 B-F3, PubMed).")
add("cat-32", M,
    "| Liver eQTL meta-analysis | Broadaway et al. (2024) [ref] |",
    "| Liver eQTL meta-analysis | Broadaway et al., 2024 |", ["B-F3"], "Ref 10.")
add("cat-33", M,
    "| Liver caQTL | Currin et al. (2025) [ref] |",
    "| Liver caQTL | Currin et al., 2025 |", ["B-F3"], "Ref 9.")
add("cat-34", M,
    "| Liver ABC enhancers | [ref: liver ABC-enhancer source] |",
    "| Liver ABC enhancers | Nasser et al., 2021 |", ["B4"], "ABC maps of Nasser et al. 2021 (PMID 33828297).")
add("cat-35", M,
    "| Histology-anchored continuum signature | Kamzolas et al. (2026) [ref] |",
    "| Histology-anchored continuum signature | Kamzolas et al., 2026 |", ["B-F3"], "Ref 1.")
add("cat-36", M,
    "Weiss et al. (2026), Supplementary Table 3 [ref]",
    "Weiss et al., 2026, Supplementary Table 3", ["B-F3"], "Ref 20.")
add("cat-37", M,
    "| Published MASLD gene lists | Kumak et al. (2026) [ref] |",
    "| Published MASLD gene lists | Kumak et al., 2026 |", ["B-F3"], "Ref 14.")

# ---------------- B-F2: InstaPrism ----------------
add("cat-38", M,
    "| BayesPrism and InstaPrism | [ref: BayesPrism]; [ref: InstaPrism] | InstaPrism 0.1.6; BayesPrism [version to be supplied] |",
    "| InstaPrism (inference under the BayesPrism model) | [ref: BayesPrism]; [ref: InstaPrism] | InstaPrism 0.1.6 |",
    ["B-F2"], "Every fraction object that feeds a reported analysis is InstaPrism class; no BayesPrism run feeds one. The 0.1.6 version is not confirmed by the environment (held).",
    [])
add("cat-39", M,
    "two fraction inputs (BayesPrism fractions estimated from the fragment counts used here and from an earlier count matrix of the same cohorts)",
    "two fraction inputs (InstaPrism estimates under the BayesPrism model, from the fragment counts used here and from an earlier count matrix of the same cohorts)",
    ["B-F2"], "class(readRDS) = InstaPrism for corrected_774 and earlier results/ objects (GSE126848, GSE162694).")
add("cat-40", M,
    "from the 16 estimated BayesPrism lineage fractions [ref: BayesPrism].",
    "from the 16 lineage fractions estimated with InstaPrism under the BayesPrism model [ref: BayesPrism] [ref: InstaPrism].",
    ["B-F2"], "The continuum composition inputs are Analysis/Deconvolution/results/{GSE162694,GSE213621}/*_bayesprism_result.rds, both InstaPrism class (checked by cat).")
add("cat-41", M,
    "(fragment-count and earlier-matrix BayesPrism fractions)",
    "(fragment-count and earlier-matrix InstaPrism fractions)",
    ["B-F2"], "Same estimator naming in the Fig 4E Methods.")
add("cat-42", D + "fig4.md",
    "three fibrosis calls remain with BayesPrism cell-type fractions",
    "three fibrosis calls remain with InstaPrism cell-type fractions",
    ["B-F2"], "Same estimator naming in the Fig 4E text.")

# ---------------- B-F3: lead contact, duplicate markers ----------------
add("cat-43", M,
    "Requests for further information should be directed to the lead contact, [name and email to be supplied].",
    "Requests for further information should be directed to the lead contact, Neville E. Sanjana (nsanjana@nygenome.org).",
    ["B-F3"], "Matches title.md:12.")
add("cat-44", M,
    "GSE296875, GSE192741, GSE287826 and HRA007511 (HMSMA).",
    "GSE296875, GSE192741, the Vu et al. (2026) Visium deposit (UQ eSpace, doi:10.48610/e95155f), the Yakubovsky et al. (2026) Visium deposit (Zenodo, doi:10.5281/zenodo.17735506 and doi:10.5281/zenodo.17735587) and HRA007511 (HMSMA; Li et al., 2025). The GeoMx series GSE287826 was obtained but not analyzed.",
    ["B-F1"], "Data availability must name the analyzed spatial sources and not list GSE287826 as an input.")
add("cat-45", M, "(McCarthy and Smyth, 2009) [ref: McCarthy and Smyth, 2009]", "(McCarthy and Smyth, 2009)", ["B-F3"], "Duplicate marker; ref 37.")
add("cat-46", M, " F0 [ref: Hoang et al., 2019]", " F0", ["B-F3"], "The sentence already names Hoang et al. (2019).")
add("cat-47", M, "(Gerard, 2020) [ref: Gerard, BMC Bioinformatics 2020]", "(Gerard, 2020)", ["B-F3"], "Ref 39.")
add("cat-48", M, "Kamzolas et al. [ref: Kamzolas et al., 2026; doi:10.1038/s42255-026-01543-7]", "Kamzolas et al. (2026)", ["B-F3"], "Ref 1.")
add("cat-49", M, "those of Kamzolas et al. [ref: Kamzolas et al., 2026]", "those of Kamzolas et al. (2026)", ["B-F3"], "Ref 1.")
add("cat-50", M, "to the Gene Catalog [ref: Kumak et al., 2026]", "to the Gene Catalog", ["B-F3"], "The sentence already names Kumak et al. (2026).")
add("cat-51", M, "Methods of Filliol et al. [ref: Filliol et al.; PMC9949942]", "Methods of Filliol et al. (2022)", ["B-F3"], "Ref 15.")
add("cat-52", M, "Guilliams et al. Table S6 [ref: Guilliams et al., Table S6]", "Guilliams et al. (2022), Table S6", ["B-F3"], "Ref 16.")
add("cat-53", M, "(REML; `nlme::gls` with `varIdent`) [ref: nlme]", "(REML; `nlme::gls` with `varIdent`; Pinheiro et al.)", ["B-F3"], "Ref 38.")
add("cat-54", M, "GSE192741 and Vu et al. [ref: Vu et al.] test", "GSE192741 and Vu et al. (2026) test", ["B-F3"], "Ref 17.")
add("cat-55", M, "tissue-context analysis [ref: Yakubovsky et al., 2026].", "tissue-context analysis (Yakubovsky et al., 2026).", ["B-F3"], "Ref 19.")
add("cat-56", M, "(M4 to M8; unit the donor) [ref: Yakubovsky et al., 2026]", "(M4 to M8; unit the donor)", ["B-F3"], "The sentence already names Yakubovsky et al. (2026).")
add("cat-57", M, "and was not refitted [ref: Weiss et al., 2026]", "and was not refitted", ["B-F3"], "The sentence already names Weiss et al. (2026).")
add("cat-58", M, "(Boesch et al., 2026) [ref: Boesch et al., 2026]", "(Boesch et al., 2026)", ["B-F3"], "Duplicate marker.")
add("cat-59", D + "fig4.md", "[ref: Kumak et al., 2026, bioRxiv, doi:10.64898/2026.08.26.747338]", "(Kumak et al., 2026)", ["B-F3"], "Ref 14.")
add("cat-60", D + "fig3.md", "[ref: Kumak et al., 2026, bioRxiv, doi:10.64898/2026.08.26.747338]", "(Kumak et al., 2026)", ["B-F3"], "Ref 14.")
add("cat-61", D + "fig5_discussion.md", "(Kumak et al., 2026 [ref: bioRxiv\ndoi:10.64898/2026.08.26.747338])", "(Kumak et al., 2026)", ["B-F3"], "Ref 14.")
add("cat-62", D + "fig5_discussion.md", "(Kamzolas et al., 2026 [ref: doi:10.1038/s42255-026-01543-7])", "(Kamzolas et al., 2026)", ["B-F3"], "Ref 1.")
add("cat-63", D + "fig5_discussion.md", "(Gobeil et al., 2026b [ref: medRxiv doi:10.64898/2026.06.04.26354903])", "(Gobeil et al., 2026b)", ["B-F3"], "Ref 13.")
add("cat-64", D + "fig5_discussion.md", "[ref: Oakes and Pachter, doi:10.1186/s12863-026-01440-8]", "(Oakes and Pachter, 2026)", ["B-F3"], "Ref 27.")
add("cat-65", L, "Zhu et al. [ref: Zhu et al., liver MPRA of disease-associated variants])", "Zhu et al., 2026)", ["B-F3"], "Ref 40.")
add("cat-66", L, "Elison et al. [ref: Elison et al., liver single-cell chromatin and caQTL study])", "Elison et al., 2025)", ["B-F3"], "Ref 25.")
add("cat-67", L, "Vu et al. [ref: Vu et al., liver spatial transcriptomics]", "Vu et al. (2026)", ["B-F3"], "Ref 17.")
add("cat-68", D + "fig2.md", "(0.706, European NAFLD GWAS [ref: 2021 European NAFLD GWAS, PMID 34841290])", "(0.706, European NAFLD GWAS; Ghodsian et al., 2021)", ["B-F3"], "Cites ref 11 by name where its study is used.")
add("cat-69", D + "fig2.md", "[ref: UK Biobank liver-enzyme GWAS]", "(Sinnott-Armstrong et al., 2021)", ["B-F3"], "Ref 35 (METHODS.md:34).")
add("cat-70", D + "fig2.md", "liver ABC-enhancer [ref: liver ABC enhancer-gene maps], and", "liver ABC-enhancer (Nasser et al., 2021), and", ["B4"], "Adds the ABC source.")

# ---------------- GSE287826 reconciliation ----------------
add("cat-71", M,
    "GSE287826 was not analyzed, because no authoritative donor key was public.",
    "GSE287826, a NanoString GeoMx series, was not analyzed (see Assay observability of programs).",
    ["B-F1"], "One stated reason, given once at :427 (critic C14 missed[0]).")
add("cat-72", M,
    "GeoMx was excluded from inference because no gene axis was loaded and its direction and normalization were unresolved. GSE287826 (no gene axis loaded; no AOI-to-donor key) tests no complete program in the current processing, so its observability is not yet determined.",
    "GSE287826, a NanoString GeoMx series of hepatocyte areas of interest from 4 healthy and 11 MASH livers, was obtained but not analyzed. Its gene axis was not loaded, so its normalization, direction and program observability were not determined, and no public key joins its areas of interest to donors.",
    ["B-F1"], "Merges the two descriptions of the same GeoMx series into one.",
    ["Analysis/Spatial/data/public_expansion_2026/gse287826/GSE287826_family.soft.gz: Series_overall_design '4 healthy and 11 MASH patients ... NanoString GeoMx'"])

# ---------------- B-F4: years, Gobeil a/b, LGALS1 ----------------
add("cat-73", D + "references.md", "*J Clin Endocrinol Metab* **111**, 1980–1993 (2026).", "*J Clin Endocrinol Metab* **111**, 1980–1993 (2026a).", ["B-F4"], "Matches the 2026a citation (01_intro.md:7).")
add("cat-74", D + "references.md", "https://doi.org/10.64898/2026.06.04.26354903 (2026).", "https://doi.org/10.64898/2026.06.04.26354903 (2026b).", ["B-F4"], "Matches the 2026b citations.")
add("cat-75", D + "fig5.md", "Vu et al., 2025)", "Vu et al., 2026)", ["B-F4"], "Ref 17 is JHEP Rep 8, issue 2 (2026); PubMed 41541503.")
add("cat-76", D + "fig5_discussion.md",
    "reach that set (best ABF PP.H4 0.47; Catalog class unresolved) and LGALS1 in the established-state class.",
    "reach that set (best ABF PP.H4 0.47; Catalog class unresolved), and LGALS1 is in the established-state class.",
    ["B-F4"], "Adds the missing verb; LGALS1 is established_state_associated in r9d (checked by cat).")

# ---------------- references: remove uncited 41, add verified entries ----------------
NEWREFS = (
    "41. Li, Z. et al. Spatially resolved multi-omics of human metabolic dysfunction-associated steatotic liver disease. *Nat Genet* **57**, 3112–3125 (2025). https://doi.org/10.1038/s41588-025-02407-8\n\n"
    "42. Hong, S. E. et al. Single-cell eQTL analysis identifies genetic variation underlying metabolic dysfunction-associated steatohepatitis. *Nat Genet* **57**, 1638–1648 (2025). https://doi.org/10.1038/s41588-025-02237-8\n\n"
    "43. Kim, H. Y. et al. Multi-modal analysis of human hepatic stellate cells identifies novel therapeutic targets for metabolic dysfunction-associated steatotic liver disease. *J Hepatol* **82**, 882–897 (2025). https://doi.org/10.1016/j.jhep.2024.10.044\n\n"
    "44. Nasser, J. et al. Genome-wide enhancer maps link risk variants to disease genes. *Nature* **593**, 238–243 (2021). https://doi.org/10.1038/s41586-021-03446-x\n\n"
    "45. Seagle, H. M. et al. Genomics-informed drug-repurposing strategy identifies two therapeutic targets for preventing liver disease associated with metabolic dysfunction. *Am J Hum Genet* **112**, 1778–1791 (2025). https://doi.org/10.1016/j.ajhg.2025.06.014\n\n"
    "46. Kozumi, K. et al. Transcriptomics identify thrombospondin-2 as a biomarker for NASH and advanced liver fibrosis. *Hepatology* **74**, 2452–2466 (2021). https://doi.org/10.1002/hep.31995\n\n"
    "47. Lefebvre, P. et al. Interspecies NASH disease activity whole-genome profiling identifies a fibrogenic role of PPARα-regulated dermatopontin. *JCI Insight* **2**, e92264 (2017). https://doi.org/10.1172/jci.insight.92264"
)
add("cat-77", D + "references.md",
    "41. Myojin, Y. et al. Genome-wide CRISPR screen identifies ACSL3 as a regulator of lipotoxicity and progression of MASLD. *Hepatol Commun* **10**, e0884 (2026). https://doi.org/10.1097/HC9.0000000000000884",
    NEWREFS, ["B-F3", "B4"],
    "Removes Myojin (cited nowhere; its data enter only the terminal Plan 40-44 stress tests, RESULTS.md:5572) and adds PubMed-verified source studies now cited by name.",
    ["PubMed PMIDs 41286103, 40562914, 39522884, 33828297, 40780050, 34105780, 28679947 (volume/pages from PubMed metadata)"])

# ---------------- A-NOV precedent (needed so Seagle/Hong/Zhu are cited) ----------------
add("cat-78", D + "01_intro.md",
    "None of these studies keeps fine-mapped multi-ancestry colocalization, multicohort disease-state expression",
    "None of these studies keeps colocalization of fine-mapped multi-ancestry GWAS with a European liver eQTL panel, multicohort disease-state expression",
    ["A-NOV"], "Only the GWAS side spans ancestries (METHODS.md:248, :252).")
add("cat-79", D + "fig5_discussion.md",
    "expression and coding evidence to nominate genes (Sveinbjornsson et al., 2022).",
    "expression and coding evidence to nominate genes (Sveinbjornsson et al., 2022).\nMulti-ancestry GWAS have been combined with genetically predicted liver\nexpression to nominate 212 candidate MASLD genes (Seagle et al., 2025).\nSingle-cell eQTL mapping in liver biopsies from 25 people with MASLD and 23\ncontrols found cell-state-dependent eQTLs enriched for MASLD heritability (Hong\net al., 2025). Liver chromatin accessibility, a reporter assay, liver eQTL,\nchromatin looping and CRISPRi have been combined to assign target genes to\nMASLD risk variants (Zhu et al., 2026).",
    ["A-NOV"], "Credits the nearest precedents; Zhu is ref 40 but was cited only in a legend.",
    ["PubMed 40780050 abstract (212 genes)", "PubMed 40562914 abstract (25 MASLD, 23 controls)", "PubMed 42230773 abstract"])

# ---------------- checks ----------------
ok = True
texts = {}
for e in E:
    p = ROOT / e["file"]
    t = texts.setdefault(e["file"], p.read_text())
    n = t.count(e["old"])
    if n != 1:
        ok = False
        print(f"COUNT {n}: {e['id']} {e['file']} :: {e['old'][:80]!r}")
# overlap check
by_file = {}
for e in E:
    t = texts[e["file"]]
    s = t.find(e["old"]); by_file.setdefault(e["file"], []).append((s, s + len(e["old"]), e["id"]))
for f, spans in by_file.items():
    spans.sort()
    for a, b in zip(spans, spans[1:]):
        if b[0] < a[1]:
            ok = False; print("OVERLAP", f, a[2], b[2])
# simple style checks on new text
for e in E:
    if "—" in e["new"]:
        ok = False; print("EM DASH", e["id"])
    if re.search(r"\bTREAT\b|\bcausal gene\b|progression", e["new"]):
        print("WORD CHECK", e["id"], re.findall(r"\bTREAT\b|\bcausal gene\b|progression", e["new"]))
print("edits", len(E), "ok", ok)
if not ok:
    sys.exit(1)

# apply in memory and list remaining [ref placeholders
applied = dict(texts)
for e in E:
    applied[e["file"]] = applied[e["file"]].replace(e["old"], e["new"], 1)
src = ["docs/manuscript/draft/" + x for x in ["abstract.md", "title.md", "01_intro.md", "fig1.md", "fig2.md", "fig3.md", "fig4.md", "fig5.md", "fig5_discussion.md", "fig7.md", "supplementary_notes.md"]] + [M, L]
rows = []
for f in src:
    t = applied.get(f) or (ROOT / f).read_text()
    for m in re.finditer(r"\[(ref[^\]]*|[^\]]*to be (?:supplied|assigned|confirmed)[^\]]*)\]", t):
        ln = t.count("\n", 0, m.start()) + 1
        ctx = t[max(0, m.start() - 60):m.start()].replace("\n", " ").replace("\t", " ")
        rows.append((f, ln, m.group(0).replace("\n", " "), ctx))
out = RUN / "staged/ref_placeholders.tsv"
with out.open("w") as fh:
    fh.write("file\tline_after_cat_edits\tplaceholder\tpreceding_text\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print("placeholders", len(rows), "ref-type", sum(r[2].startswith("[ref") for r in rows))
(RUN / "staged/edits_cat.json").write_text(json.dumps(E, indent=1, ensure_ascii=False))
print("wrote", RUN / "staged/edits_cat.json")
