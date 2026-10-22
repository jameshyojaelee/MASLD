# Builds the proposed text edits for the Methods workstream. old_text is read from the
# files by line range so it matches exactly; nothing in docs/ is modified.
import json, sys
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
OUT = sys.argv[1]
P11_04 = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None


def lines(f, a, b):
    with open(R + f) as h:
        L = h.read().split("\n")
    return "\n".join(L[a - 1:b])


def ed(fid, f, a, b, new, claim=False, sub=None):
    old = lines(f, a, b)
    if sub is not None:
        assert sub in old, (fid, f, a, b, sub)
        old = sub
    return dict(fid=fid, file=f, line=f"{a}-{b}" if a != b else str(a), old_text=old, new_text=new, claim_level=claim)


E = []
M = "docs/manuscript/METHODS.md"
# ---------------- R10-02
E.append(ed("R10-02", "docs/manuscript/draft/fig5_discussion.md", 76, 78,
"a target-validity claim. Its failed matched clinical benchmark and five public\n"
"perturbation stress tests that did not meet their prespecified criteria (Methods)\n"
"remain falsification boundaries. The publication-facing portal must default to the"))
E.append(ed("R10-02", M, 810, 810,
lines(M, 810, 810) + "\n\n"
"## Supplementary public perturbation stress tests\n\n"
"Five prespecified analyses asked whether Resource classes or programs predicted\n"
"responses in public perturbation or follow-up models. (i) Evidence classes and\n"
"palmitate-lipotoxicity CRISPR survival in the Myojin hepatocyte screen: the\n"
"proposed Resource-level enrichment was not confirmed; survival is not\n"
"interchangeable with lipid storage. (ii) Induction, localization and reversal of\n"
"the two frozen hepatocyte programs by disease cues in precision-cut liver slices\n"
"and organoids: the multi-assay criteria did not pass. (iii) A double dissociation\n"
"in which genetic-only genes are functional and state-only genes are spatially\n"
"organized: the two-arm criteria did not pass. (iv) Reversal of the 117-program\n"
"axis with histologic improvement: the risk-to-state criteria did not pass.\n"
"(v) A lineage-axis topology of fat-versus-profibrotic multicellular response:\n"
"no source-independent result. None of these analyses defines an evidence class,\n"
"threshold, hero gene or figure; they bound what acute public models can test."))
# ---------------- R11-05 (gloss of source_dependent for ATAC)
E.append(ed("R11-05", "docs/RESULTS.md", 2698, 2698,
"| ATAC `source_dependent` peaks | 1,844 / 39,914 | BH q < 0.05 in one cohort only; 1,505 of 1,844 have the same logFC sign in the other cohort, so most reflect unequal power (6/6 versus 5/9 donors), not disagreement |"))
E.append(ed("R11-05", "docs/RESULTS.md", 2906, 2907,
"Across the 39,914 jointly testable peaks, evidence states were 4 `supported`,\n"
"1,844 `source_dependent` (q < 0.05 in one cohort only; same sign in the other\n"
"cohort for 1,505), and 38,066 `indeterminate`, with no `discordant`"))
# ---------------- R11-06 (assay limitation paragraph) + R11-05 gloss inside it
E.append(ed("R11-06", "docs/PAPER.md", 312, 318,
"Observability limits are an important Resource result, and they have three\n"
"different causes. Panel content limits CosMx, which adequately observes 12 of\n"
"117 programs. Unresolved data processing limits GeoMx (direction and\n"
"normalization unresolved) and GSE287826 (no AOI-to-donor key); each currently\n"
"tests 0 of 117 programs, which is not a property of the assay. Donor number\n"
"limits single-cell ATAC (5 normal and 9 MASH donors in GSE244832, 6 and 6 in\n"
"GSE281367): 38,066 of 39,914 jointly testable peaks are `indeterminate`, and\n"
"1,844 are `source_dependent` (q < 0.05 in one cohort only; 1,505 of them have\n"
"the same sign in the other cohort). Against these, the four Visium/ST datasets\n"
"each observe 114 of 117 programs and the liver proteome meets the coverage rule\n"
"for 26 of 117. The current Catalog has zero strict `tested_negative` calls.\n"
"These are not failed validations. They show why modality-vote counting converts\n"
"absent coverage, unresolved donor identity, and assay mismatch into false\n"
"confidence.", claim=False))
E.append(ed("R11-06", M, 1079, 1085,
"Complete per-assay coverage is in the Supplement. In summary: CosMx adequately\n"
"observes 12 of 117 programs because of its targeted panel; GeoMx (direction and\n"
"normalization unresolved) and GSE287826 (no AOI-to-donor key) test no complete\n"
"program in the current data processing; the four Visium/ST datasets each\n"
"observe 114 of 117 programs; and the liver proteome meets the coverage rule for\n"
"26 of 117. In single-cell ATAC, with 5/9 and 6/6 normal/MASH donors per cohort,\n"
"38,066 of 39,914 jointly testable peaks are `indeterminate` and 1,844 are\n"
"`source_dependent` (q < 0.05 in one cohort only; 1,505 share sign with the other\n"
"cohort); and the Gene Catalog contains zero strict `tested_negative` calls.\n"
"These are limits on what the assembled assays and data can currently resolve,\n"
"and they are the reason modality-vote counting is not used anywhere in this work."))
E.append(ed("R11-06", "docs/RESULTS.md", 2694, 2695,
"| GeoMx complete-program testability | 0 / 117 programs | `untestable`: 0 program genes on the loaded gene axis; direction and normalization unresolved (data-processing state, not an assay limit) |\n"
"| GSE287826 complete-program testability | 0 / 117 programs | `untestable`: 0 program genes on the loaded gene axis and no authoritative AOI-to-donor key (data-processing state, not an assay limit) |\n"
"| Visium/ST complete-program observability (GSE192741, Vu, HMSMA, Yakubovsky) | 114 / 117 programs each | 3 partial in each; 0 untestable |"))
# ---------------- R11-08
E.append(ed("R11-08", "docs/STATUS.md", 1142, 1145,
"- ATAC v3 separates static observability from MASH-associated accessibility.\n"
"  Four stellate-lineage peaks reach q < 0.05 in both cohorts with the same\n"
"  sign, against 1.33 expected by chance from the per-cohort counts (820 and 15\n"
"  of 9,264; hypergeometric P(X >= 4) = 0.038, one uncorrected test treating\n"
"  peaks as independent). They are candidates, not chromatin findings, and not\n"
"  causal-gene or causal-cell-type proof; three of the four promoters are\n"
"  unnamed ENSG loci. The GSE281367 T/NK label remains a\n"
"  combined, source-dependent observability context."))
E.append(ed("R11-08", "docs/RESULTS.md", 2908, 2912,
"peak. All four supported peaks were stellate-lineage peaks, retained direction\n"
"under both leave-one-donor-out analyses, and did not reverse under the\n"
"GSE244832 batch or broad MASLD sensitivities. Their two-cohort overlap is\n"
"4 observed against 1.33 expected (820 and 15 significant of 9,264;\n"
"hypergeometric P(X >= 4) = 0.038, uncorrected, peaks treated as independent).\n"
"This satisfies the prespecified\n"
"candidate dynamic-main eligibility rule, but it is not yet an adopted main-panel\n"
"claim and does not prove a causal gene or causal cell type."))
# ---------------- R10-04
E.append(ed("R10-04", M, 94, 94,
"Historical technical inclusion used the `pass_technical` flag, defined by PCA and library-size quality control. Samples failing only the expression-based sex check remained eligible for the reported analyses, in which sex enters only as a model covariate.",
sub=lines(M, 94, 94)))
# ---------------- R09-04: drop the unadopted allele-candidate methods
E.append(ed("R09-04", M, 953, 989, ""))
# ---------------- R16-04 (STATUS H17 row)
E.append(ed("R16-04", "docs/STATUS.md", 315, 315,
"| H17 | Null (uninformative in vivo arm) | Source-corrected refit completed 2026-10-14 (job 22237956): do iPSC and in vivo hepatocytes share a lipid-loading response? Unit: participant for iPSC-hepatocytes (193, paired); donor for in vivo hepatocytes (55: 9 control, 46 steatotic; GSE244832 and GSE202379); 14,975 genes. Same sign among 7,628 iPSC q < 0.05 genes: 0.447 (donor-label permutation 95th percentile 0.603, lower-tail p 0.21, 1,000 permutations). Spearman rho -0.093 (donor bootstrap 95% CI -0.177 to 0.036). Significant in both: 3 against 3.06 expected (p 0.67); core 1 gene (MCTP2). All five prespecified support criteria failed; iPSC positive control passed. The in vivo arm calls 6 genes at q < 0.05 and misses CIDEC (q 0.71), PLIN2 (q 0.9985) and PLIN4 (q 0.63), so it cannot detect canonical lipid-droplet genes. No opposite-direction claim: the gene-shuffle p 1e-4 treats correlated genes as independent and the steatosis-only sensitivity reverses sign. | Bulk signature from the 774-unit refit (adopted 2026-10-11); in vivo stage labels from the source-corrected roster | `H17-source-stage-corrected-20261014_22237956/decision.tsv` |"))
# ---------------- R16-07 authority staleness
E.append(ed("R16-07", "docs/ROADMAP.md", 11, 15,
"installed. The source-diagnosis H17 and N3 refits completed on 2026-10-14\n"
"(H17: no shared lipid-loading response, in vivo arm uninformative; N3: no\n"
"qualified corrected calls). Remaining publication work includes final\n"
"publication layout and the overall checklist below; unavailable\n"
"donor/histology joins remain explicit scientific limits."))
E.append(ed("R16-07", "docs/ROADMAP.md", 136, 137,
"The corrected COLOC rerun adopted on 2026-08-17 is superseded by COLOC r2, adopted\n"
"on 2026-10-14; its signal eligibility and trait drivers were corrected (see the",
))
E.append(ed("R16-07", "docs/ROADMAP.md", 151, 152,
"wording corrections are complete as candidates; the 774-unit bulk refit was\n"
"adopted on 2026-10-11 and COLOC r2 on 2026-10-14."))
E.append(ed("R16-07", "docs/STATUS.md", 458, 458,
"of the refit reproduction, 24,115 rows), and on COLOC r2, adopted 2026-10-14 after r8 was built. The user approved the"))
E.append(ed("R16-07", "docs/STATUS.md", 828, 828,
lines("docs/STATUS.md", 828, 828).replace("from the adopted 2026-08-17 COLOC release", "from COLOC r2 (adopted 2026-10-14)")))
E.append(ed("R16-07", "docs/STATUS.md", 976, 976,
"   DNA tables must now be rebuilt from COLOC r2 (adopted 2026-10-14)."))
E.append(ed("R16-07", "docs/STATUS.md", 999, 999,
"refit with COLOC r2 (adopted 2026-10-14) `susie_coloc_r2_20260924T053244Z`"))
E.append(ed("R16-07", "docs/STATUS.md", 1094, 1097,
"Exact genetic counts and overlap percentages are now authorized **from COLOC r2,\n"
"adopted 2026-10-14**, as ledgered in [RESULTS.md](RESULTS.md); the July and\n"
"2026-08-17 values are superseded. Gene examples and Gene Catalog classes use\n"
"their r2 rebuilds."))
E.append(ed("R16-07", "docs/RESULTS.md", 2051, 2052,
"**Do not upgrade the depletion to a finding.** On the 774 units with COLOC r2\n"
"(adopted 2026-10-14) the overlap test gives OR 0.76, nominal p = 0.20 in the"))
E.append(ed("R16-07", "docs/RESULTS.md", 2221, 2221,
"with COLOC r2, adopted 2026-10-14, 92.2% expected against 93.9% observed for 426, 1,155"))
E.append(ed("R16-07", "docs/README.md", 125, 125,
"summaries, not inferential or temporal bins. The exact GNMT direct-NAFLD ABF pair (PP.H4 0.7058, 2021_34841290_NAFLD_EUR) is in COLOC r2; its local haplotype remains unresolved."))
E.append(ed("R16-07", "docs/PAPER.md", 612, 613,
"- A preliminary ambient-RNA audit raises a lineage-assignment concern for the\n"
"  Stromal ECM and Ductular injury programs. An aligned rerun is independently\n"
"  validated but not adopted. Their"))
# ---------------- R15-07 data citation
E.append(ed("R15-07", M, 558, 558,
"The primary genetic portfolio contained 35 GWAS strata (Supplementary Table GWAS-sources): 12 from the Million Veteran Program Release 4 genetically inferred ancestry summary statistics (open-access dbGaP phs002453.v1.p1; NAFLD, ALT and AST in AFR, AMR, EAS and EUR), 9 from published GWAS (PMIDs 31311600, 32298765, 34128465, 34841290, 34957434, 36402844 and 36280732, the last contributing deCODE, Intermountain and UK Biobank strata), 2 from FinnGen R12 (NAFLD, NASH), 3 UK Biobank European liver-enzyme strata, 3 Biobank Japan liver-enzyme strata and 6 Pan-UK Biobank AFR and CSA liver-enzyme strata. Strata are `placement = main` and `tier` 1 or 2 in `GWAS/finemapping/config/gwas_trait_tier.tsv`."))
E.append(ed("R15-07", M, 1024, 1024,
"This study integrates publicly available datasets, cited with their accessions in the relevant Methods subsections, Supplementary Table GWAS-sources, source registry, and reference list. Single-cell, ATAC and spatial inputs include GSE192742 (Guilliams liver atlas), GSE244832 and GSE281367 (snATAC), GSE296875, GSE192741, GSE287826 and HRA007511 (HMSMA). The final Resource release will include source checksums, access dates, testing-universe manifests, a data dictionary, and the source-dependency graph used by the Gene Catalog."))
# ---------------- R15-08 code availability
E.append(ed("R15-08", M, 1028, 1040,
"The current numbers are generated by the following files.\n\n"
"- Bulk 774-unit refit: `scripts/manuscript/resource_f_five/refit_corrected_controls.R`\n"
"- COLOC r2: `GWAS/finemapping/src/06_susie_coloc.R`, `GWAS/finemapping/src/07_combine_susie_coloc.R`\n"
"- Figure 2 inference: `scripts/figures/fig2_trait_directness_stats.py`, `scripts/figures/fig2_multiplicity_correction.py`\n"
"- Gene Catalog r8: `scripts/portal/build_evidence_passport_bundle.py`, `scripts/portal/validate_evidence_passports.py`\n\n"
"Historical components (supplementary only): `scripts/manuscript/build_evidence_class_release.R`,\n"
"`scripts/manuscript/build_strict_drug_benchmark.py`, `Analysis/Spatial/scripts/15l_gsmap_pcc_narrowing.py`,\n"
"`scripts/manuscript/build_gsmap_class_validation.py`, `scripts/manuscript/build_coding_architecture_release.R`,\n"
"`scripts/portal/generate_evidence_passports.py`. Their `passport` filenames are retained because they are\n"
"hash-bound to the recorded v1 candidate; the manuscript and portal call the product the MASLD Gene Catalog."))
if P11_04:
    E += P11_04
json.dump(E, open(OUT, "w"), indent=1, ensure_ascii=False)
print(len(E), "edits")
