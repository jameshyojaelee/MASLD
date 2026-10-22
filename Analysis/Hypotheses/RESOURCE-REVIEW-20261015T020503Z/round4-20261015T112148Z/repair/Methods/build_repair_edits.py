"""Build the repair-round text edits for the Methods workstream and check them against the current files.

Numbers for R11-04 and R11-08 are read from the repair summary table written by combine_repair.R.
Every old_text must occur exactly once in the current file; new_text must have no em dash and no
evidence tag. Output: edits_repair.json in this directory.
"""
import csv
import json
import os
import re

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HERE = os.path.join(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/Methods")
FIG = os.path.join(ROOT, "figures/candidates/review-r3-fixes-20261015T112148Z/Methods/repair")
FIRST = os.path.join(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair_inputs/Methods.json")

with open(os.path.join(FIG, "R11-04_R11-08_donor_relabel_summary_repair.tsv")) as fh:
    S = {r["lineage"]: r for r in csv.DictReader(fh, delimiter="\t")}
hep, st = S["hepatocyte"], S["stellate"]


def f3(x):
    return f"{float(x):.3f}"


def f4(x):
    return f"{float(x):.4f}"


N = {
    "hep_obs": int(hep["observed_same_sign"]), "hep_n": int(hep["n_ge_observed_same_sign"]), "hep_p": f4(hep["p_same_sign"]),
    "hep_sd": f"{float(hep['null_sd_same_sign']):.0f}", "hep_mean": f"{float(hep['null_mean_same_sign']):.0f}",
    "st_n": int(st["n_ge_observed_same_sign"]), "st_p": f4(st["p_same_sign"]),
    "q3_hep": f3(hep["q_BH_three_rows_release_family"]), "q3_st": f3(st["q_BH_three_rows_release_family"]),
    "q2_hep": f4(hep["q_BH_two_testable_rows_post_hoc"]), "q2_st": f4(st["q_BH_two_testable_rows_post_hoc"]),
    "sup_n": int(st["n_ge_observed_supported"]), "sup_p": f4(st["p_supported"]),
    "sup_mean": f"{float(st['null_mean_supported']):.2f}",
}
assert N["hep_obs"] == 849
print(N)

first = {e["fid"]: e for e in json.load(open(FIRST))["entries"]}


def te(file, line, old, new):
    return {"file": file, "line": line, "old_text": old, "new_text": new, "claim_level": False}


def keep(fid, i):
    t = dict(first[fid]["text_edits"][i])
    return t


E = {}

# ---------------- R11-04 ----------------
old_tab = first["R11-04"]["text_edits"][0]["old_text"]
q3 = N["q3_st"]
new_tab = (
    "| Lineage | Donors, GSE244832 normal/MASH | Donors, GSE281367 normal/MASH | Joint peaks | Significant, GSE244832/GSE281367 | Strict replicated | Directional replication: same sign / GSE281367 discoveries; donor-relabeling p; BH q over three lineage rows | Main dynamic criterion |\n"
    "|---|---:|---:|---:|---:|---:|---:|---|\n"
    f"| Hepatocyte | 5 / 9 | 6 / 6 | 30,650 | 1 / 1,016 | 0 | 849 / 1,016; p = {N['hep_p']}; q = {N['q3_hep']} | Did not pass the sensitivity-stable replicated-peak criterion |\n"
    f"| Stellate | 5 / 9 | 6 / 6 | 9,264 | 820 / 15 | 4 | 15 / 15; p = {N['st_p']}; q = {N['q3_st']} | Did not pass: directional q = {N['q3_st']} is above 0.05 |\n"
    "| Macrophage | 4 / 8 | 6 / 6 | 0 | 0 / 0 | 0 | not testable; p = 1 | Untestable after default filtering |\n"
    "\n"
    "Directional replication is calibrated against all 2,002 relabelings of the 14\n"
    "GSE244832 donors (9 MASH), refitting the release model each time\n"
    "(condition-blind `filterByExpr`, TMM, `voomWithQualityWeights`, robust eBayes);\n"
    "the GSE281367 fits are fixed. The identity relabeling reproduces the saved\n"
    "GSE244832 log fold changes (maximum absolute difference 5e-15). Exact p is the\n"
    "fraction of relabelings, identity included, with at least the observed same-sign\n"
    f"count: hepatocyte {N['hep_n']} / 2,002, stellate {N['st_n']} / 2,002. BH adjustment uses the\n"
    "prespecified family of three lineage rows in `04_run_da.R`, with untestable\n"
    "macrophage at p = 1. Over the two testable lineages only, q would be\n"
    f"{N['q2_st']} for both; that family was not prespecified and is not used for the\n"
    "eligibility call. The earlier binomial q values (2.64×10⁻¹¹⁰ and 4.58×10⁻⁵)\n"
    "treated peaks as independent trials and are withdrawn: under donor relabeling\n"
    f"the hepatocyte same-sign count has mean {N['hep_mean']} and SD {N['hep_sd']} peaks, against SD 16\n"
    "under peak independence."
)
hip_old = "849/1,016 discovery peaks agree in direction; directional q=2.64×10⁻¹¹⁰."
hip_new = (f"849/1,016 discovery peaks agree in direction; exact GSE244832 donor-relabeling p = {N['hep_p']} "
           f"({N['hep_n']} of 2,002 relabelings), BH q = {N['q3_hep']} over the three lineage rows. The earlier "
           "binomial q (2.64×10⁻¹¹⁰) treated peaks as independent and is withdrawn.")
E["R11-04"] = [te("docs/RESULTS.md", "2900-2904", old_tab, new_tab),
               te("docs/HIGH_IMPACT_PATTERNS.md", "65", hip_old, hip_new)]

# ---------------- R11-05 ----------------
r05 = keep("R11-05", 0)
r05["new_text"] = ("| ATAC `source_dependent` peaks | 1,844 / 39,914 | BH q < 0.05 in one cohort only; 1,505 of 1,844 "
                   "have the same logFC sign in the other cohort. No peak reaches q < 0.05 in both cohorts with "
                   "opposite signs (0 `discordant`) |")
E["R11-05"] = [r05, keep("R11-05", 1)]

# ---------------- R11-06 ----------------
paper_old = first["R11-06"]["text_edits"][0]["old_text"]
paper_new = (
    "Observability limits are an important Resource result, and their causes differ.\n"
    "CosMx adequately observes 12 of 117 programs; its targeted panel covers 968\n"
    "genes. GeoMx and GSE287826 currently test 0 of 117 programs because no gene axis\n"
    "has been loaded for either (GeoMx direction and normalization are unresolved, and\n"
    "GSE287826 has no AOI-to-donor key), so their observability is not yet\n"
    "determined. The single-cell ATAC hepatocyte and stellate contrasts have 5 normal\n"
    "and 9 MASH donors in GSE244832 and 6 and 6 in GSE281367; 38,066 of 39,914\n"
    "jointly testable peaks are `indeterminate`,\n"
    "and 1,844 are `source_dependent` (q < 0.05 in one cohort only; 1,505 of them have\n"
    "the same sign in the other cohort). Each of the four Visium/ST datasets observes\n"
    "114 of 117 programs. The current Catalog has zero strict `tested_negative` calls.\n"
    "These are not failed validations. They show why modality-vote counting converts\n"
    "absent coverage, unresolved donor identity, and assay mismatch into false\n"
    "confidence."
)
meth_old = first["R11-06"]["text_edits"][1]["old_text"]
meth_new = (
    "Complete per-assay coverage is in the Supplement. In summary: CosMx adequately\n"
    "observes 12 of 117 programs (968 genes on its targeted panel); GeoMx and GSE287826\n"
    "test no complete program because no gene axis has been loaded for either (GeoMx\n"
    "direction and normalization unresolved; no GSE287826 AOI-to-donor key), so their\n"
    "observability is not yet determined; and the four Visium/ST datasets each observe\n"
    "114 of 117 programs. In single-cell ATAC, with 5/9 (GSE244832) and 6/6 (GSE281367)\n"
    "normal/MASH donors in the hepatocyte and stellate contrasts, 38,066 of 39,914\n"
    "jointly testable peaks are `indeterminate` and 1,844 are\n"
    "`source_dependent` (q < 0.05 in one cohort only; 1,505 share sign with the other\n"
    "cohort); and the Gene Catalog contains zero strict `tested_negative` calls.\n"
    "These are limits on what the assembled assays and data can currently resolve,\n"
    "and they are the reason modality-vote counting is not used anywhere in this work."
)
res_old = first["R11-06"]["text_edits"][2]["old_text"]
res_new = (
    "| GeoMx complete-program testability | 0 / 117 programs | `untestable`: no gene axis loaded; direction and normalization unresolved; observability not yet determined |\n"
    "| GSE287826 complete-program testability | 0 / 117 programs | `untestable`: no gene axis loaded and no authoritative AOI-to-donor key; observability not yet determined |\n"
    "| Visium/ST complete-program observability (GSE192741, Vu, HMSMA, Yakubovsky) | 114 / 117 programs each | 3 partial in each; 0 untestable |"
)
E["R11-06"] = [te("docs/PAPER.md", "312-318", paper_old, paper_new),
               te("docs/manuscript/METHODS.md", "1079-1085", meth_old, meth_new),
               te("docs/RESULTS.md", "2694-2695", res_old, res_new)]

# ---------------- R11-08 ----------------
st_old = first["R11-08"]["text_edits"][0]["old_text"]
st_new = (
    "- ATAC v3 separates static observability from MASH-associated accessibility.\n"
    "  Four stellate-lineage peaks reach q < 0.05 in both cohorts with the same\n"
    f"  sign. In {N['sup_n']} of 2,002 GSE244832 donor relabelings at least four peaks do so\n"
    f"  (p = {N['sup_p']}, a single test added in review; null mean {N['sup_mean']}).\n"
    "  Stellate directional replication has\n"
    f"  donor-relabeling BH q = {N['q3_st']}, so the lineage does not meet the dynamic-main\n"
    "  eligibility rule. The four peaks are candidates, not chromatin findings, and\n"
    "  not causal-gene or causal-cell-type proof; three of the four promoters are\n"
    "  unnamed ENSG loci. The GSE281367 T/NK label remains a\n"
    "  combined, source-dependent observability context."
)
rr_old = first["R11-08"]["text_edits"][1]["old_text"]
rr_new = (
    "peak. All four supported peaks were stellate-lineage peaks, retained direction\n"
    "under both leave-one-donor-out analyses, and did not reverse under the\n"
    "GSE244832 batch or broad MASLD sensitivities. Under the same 2,002 GSE244832\n"
    f"donor relabelings, {N['sup_n']} gave at least four supported peaks (p = {N['sup_p']}, a single\n"
    f"test added in review; null mean {N['sup_mean']}). The prespecified candidate dynamic-main\n"
    "eligibility rule also requires\n"
    f"directional q < 0.05; with the donor-relabeling q of {N['q3_st']} (three lineage rows),\n"
    "stellate no longer meets it. The four peaks are candidates; they are not a\n"
    "main-panel claim and do not prove a causal gene or causal cell type."
)
E["R11-08"] = [te("docs/STATUS.md", "1142-1145", st_old, st_new), te("docs/RESULTS.md", "2908-2912", rr_old, rr_new)]

# ---------------- R10-02 ----------------
r02_meth = keep("R10-02", 1)
r02_meth["new_text"] = (
    "The prespecified headline criterion required at least two approved targets in the ranking universe and convergence superiority of at least 0.03 AUROC over modality count and the best constituent channel, with a paired 95% confidence interval excluding zero. That criterion did not pass.\n"
    "\n"
    "## Supplementary public perturbation stress tests\n"
    "\n"
    "Five prespecified analyses asked whether Resource evidence classes or programs\n"
    "predicted responses in public perturbation models and paired-biopsy cohorts.\n"
    "None met its prespecified criteria. (i) Whether evidence classes or the two\n"
    "primary hepatocyte programs differ in knockout effects on three-week survival of\n"
    "HLF hepatocellular-carcinoma cells under 200 micromolar palmitate relative to\n"
    "vehicle, in the Myojin et al. genome-wide CRISPR-Cas9 knockout screen\n"
    "(Supplemental Table S1): no difference was supported (evidence-class omnibus\n"
    "permutation p = 0.537 over 7,762 genes). The screen reads cell survival, which\n"
    "is not interchangeable with lipid storage. (ii) Induction, localization and reversal of the Stromal ECM and\n"
    "Ductular injury programs in donor-paired precision-cut liver slices (GSE200418),\n"
    "human liver organoids (GSE207889, GSE253380) and paired liver biopsies\n"
    "(GSE106737): neither program met the multi-assay criteria. In the slices,\n"
    "48-hour lipid addition induced neither program, while 48 versus 24 hours of\n"
    "culture increased both. (iii) A double dissociation in which genetic-only genes\n"
    "act in cell-autonomous lipid storage and disease-state-only genes are spatially\n"
    "organized, tested in the untreated Huh7 neutral-lipid screen of the CLCC1\n"
    "study, prime-edited PNPLA3 HepG2 clones (GSE313544) and isogenic iPSC-derived\n"
    "hepatocyte-like cells (GSE158182): neither arm met its criteria, and the planned\n"
    "HMSMA Visium arm was not run because no public sample-to-clinical key exists. (iv) Reversal of the frozen 117-program state axis with histologic\n"
    "improvement in paired biopsies (GSE106737 discovery; GSE83452 and GSE48452\n"
    "replication): the nominal discovery signal did not replicate, and expression\n"
    "fingerprints showed that GSE106737 and GSE83452 contain overlapping samples\n"
    "rather than independent cohorts. (v) A lineage-topology prediction for fat,\n"
    "non-parenchymal-cell and TGF-beta cues in controlled microphysiological liver\n"
    "cultures (GSE168285) and paired biopsies from 19 participants (GSE175448): the TGF-beta\n"
    "effect and the paired-human effect were both opposite to the prediction. None\n"
    "of these analyses defines an evidence class, threshold, hero gene or figure."
)
E["R10-02"] = [keep("R10-02", 0), r02_meth]

# ---------------- R15-07 ----------------
r07 = keep("R15-07", 0)
r07["new_text"] = (
    "The primary genetic portfolio contained 35 GWAS strata (Supplementary Table GWAS-sources): 12 from the Million Veteran Program Release 4 genetically inferred ancestry summary statistics (open-access dbGaP phs002453.v1.p1; NAFLD, ALT and AST in AFR, AMR, EAS and EUR), 9 from published GWAS (PMIDs 31311600, 32298765, 34128465, 34841290, 34957434, 36402844 and 36280732, the last contributing deCODE, Intermountain and UK Biobank strata), 2 from FinnGen R12 (NAFLD, NASH), 3 UK Biobank European liver-enzyme strata (Sinnott-Armstrong et al., PMID 33462484; GWAS Catalog GCST90019492, GCST90019497 and GCST90019507), 3 Biobank Japan liver-enzyme strata (Kanai et al., PMID 29403010) and 6 Pan-UK Biobank round 2 AFR and CSA liver-enzyme strata (biomarker codes 30620, 30650 and 30730). Strata are `placement = main` and `tier` 1 or 2 in `GWAS/finemapping/config/gwas_trait_tier.tsv`."
)
E["R15-07"] = [r07, keep("R15-07", 1)]

# ---------------- R15-08 ----------------
r08 = keep("R15-08", 0)
r08["new_text"] = (
    "The bulk, COLOC, Figure 2, single-cell ATAC and Gene Catalog numbers are\n"
    "generated by the following files.\n"
    "\n"
    "- Bulk 774-unit refit: `scripts/manuscript/resource_f_five/refit_corrected_controls.R`\n"
    "- COLOC r2: `GWAS/finemapping/src/06_susie_coloc.R` (per-study SuSiE-COLOC and ABF, also used for the eligibility rerun), `GWAS/finemapping/src/perf/09_prepare_coloc_eligibility_batches.py`, `GWAS/finemapping/src/perf/10_offline_eqtl_overlap_check.R`, `GWAS/finemapping/src/perf/11_assemble_coloc_r2.py`, `GWAS/finemapping/src/perf/07_combine_rerun.R` and `GWAS/finemapping/src/perf/rebuild_tier12.py`\n"
    "- Figure 2 inference: `scripts/figures/fig2_trait_directness_stats.py`, `scripts/figures/fig2_multiplicity_correction.py`\n"
    "- Single-cell ATAC differential accessibility: `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/04_run_da.R`\n"
    "- Gene Catalog r8: `scripts/portal/prepare_passport_selection.py`, `scripts/portal/adjudicate_gen_ensembl_identity.py`, `scripts/portal/build_evidence_passport_bundle.py`\n"
    "\n"
    "Producers of the program, spatial, proteome and remaining Figure 3 to 5 numbers\n"
    "are to be added to this list before release.\n"
    "\n"
    "Historical components (supplementary only): `scripts/manuscript/build_evidence_class_release.R`,\n"
    "`scripts/manuscript/build_strict_drug_benchmark.py`, `Analysis/Spatial/scripts/15l_gsmap_pcc_narrowing.py`,\n"
    "`scripts/manuscript/build_gsmap_class_validation.py`, `scripts/manuscript/build_coding_architecture_release.R`,\n"
    "`scripts/portal/generate_evidence_passports.py`, `scripts/portal/validate_evidence_passports.py`. Their `passport`\n"
    "filenames are retained because they are hash-bound to the recorded v1 candidate; the manuscript and portal\n"
    "call the product the MASLD Gene Catalog."
)
E["R15-08"] = [r08]

# ---------------- R16-04 ----------------
r04 = keep("R16-04", 0)
r04["new_text"] = (
    "| H17 | No call (in vivo arm cannot detect canonical lipid-droplet genes) | Source-corrected refit completed 2026-10-14 (job 22237956): do iPSC and in vivo hepatocytes share a lipid-loading response? Unit: participant for iPSC-hepatocytes (193, paired); donor for in vivo hepatocytes (55: 9 control, 46 steatotic; GSE244832 and GSE202379); 14,975 genes. Same sign among 7,628 iPSC q < 0.05 genes: 0.447 (donor-label permutation 95th percentile 0.603, lower-tail p 0.21, 1,000 permutations). Spearman rho -0.093 (donor bootstrap 95% CI -0.177 to 0.036). Significant in both: 3 against 3.06 expected (p 0.67); core 1 gene (MCTP2). All five prespecified support criteria failed; iPSC positive control passed. The in vivo arm calls 7 genes at q < 0.05 (6 in the 14,975-gene shared universe) and misses CIDEC (q 0.71), PLIN2 (q 0.9985) and PLIN4 (q 0.63), so it cannot detect canonical lipid-droplet genes. No opposite-direction claim: the gene-shuffle p 1e-4 treats correlated genes as independent and the steatosis-only sensitivity reverses sign. | Bulk signature from the 774-unit refit (adopted 2026-10-11); in vivo stage labels from the 2026-10-14 source-diagnosis correction roster, the correction used for Figure 4C | `H17-source-stage-corrected-20261014_22237956/decision.tsv` |"
)
E["R16-04"] = [r04]

# ---------------- R16-07 ----------------
r7 = [dict(t) for t in first["R16-07"]["text_edits"]]
assert r7[0]["file"] == "docs/ROADMAP.md" and r7[0]["line"] == "11-15"
r7[0]["new_text"] = (
    "installed. The source-diagnosis H17 and N3 refits completed on 2026-10-14\n"
    "(H17: no shared lipid-loading response detected, and the in vivo arm could not\n"
    "detect canonical lipid-droplet genes; N3: no qualified corrected calls).\n"
    "Remaining publication work includes final publication layout and the overall\n"
    "checklist below; unavailable donor/histology joins remain explicit scientific\n"
    "limits."
)
i976 = [i for i, t in enumerate(r7) if t["file"] == "docs/STATUS.md" and t["line"] == "976"][0]
r7[i976] = te(
    "docs/STATUS.md", "974-976",
    "   DNA-only candidate are resolved. Target-transcript biotype, final examples,\n"
    "   and genetics-dependent Catalog fields are no longer waiting on compute; the\n"
    "   DNA tables must now be rebuilt from the adopted 2026-08-17 release.",
    "   DNA-only candidate are resolved. On COLOC r2 (adopted 2026-10-14), the\n"
    "   target-gene companion (424 protein-coding, 12 lncRNA), Figure 2C and the\n"
    "   exact-example source rows are installed; Gene Catalog r8 was built on r2\n"
    "   before its adoption.",
)
E["R16-07"] = r7

# ---------------- R09-04 ----------------
r0904_old = json.load(open(os.path.join(HERE, "r0904_old.json")))["old"]
E["R09-04"] = [te(
    "docs/manuscript/METHODS.md", "885-988", r0904_old,
    "A separate two-model Figure 7 development (a rank-12 RNA-to-H3K27ac model with\n"
    "local correction, its acquisition and input-robustness checks, and an\n"
    "AlphaGenome allele analysis against Currin caQTL effects) has not been adopted,\n"
    "and no result from it is reported in this paper. Its code and source tables are\n"
    "in `scripts/analysis/masld_two_models/`.\n",
)]

# ---------------- R10-04 ----------------
r104 = keep("R10-04", 0)
r104["new_text"] = (
    "Historical technical inclusion used the `pass_technical` flag, defined by PCA and library-size quality control. "
    "Samples failing only the expression-based sex check remained eligible for the reported bulk analyses, which adjust "
    "for expression-inferred sex."
)
E["R10-04"] = [r104]

# ---------------- checks ----------------
TAG = re.compile(r"\[(Inference|Unverified|Speculation|Estimate|Contains Unverified Content)")
problems = []
for fid, edits in E.items():
    for t in edits:
        s = open(os.path.join(ROOT, t["file"])).read()
        c = s.count(t["old_text"])
        line = s[: s.find(t["old_text"])].count("\n") + 1 if c else None
        t["line_found"] = line
        if c != 1:
            problems.append(f"{fid} {t['file']} old_text count {c}")
        if "—" in t["new_text"]:
            problems.append(f"{fid} em dash")
        if TAG.search(t["new_text"]):
            problems.append(f"{fid} evidence tag")
        print(fid, t["file"], t["line"], "found at", line, "count", c)
# overlapping edits within one file
by_file = {}
for fid, edits in E.items():
    for t in edits:
        s = open(os.path.join(ROOT, t["file"])).read()
        a = s.find(t["old_text"])
        by_file.setdefault(t["file"], []).append((a, a + len(t["old_text"]), fid))
for f, spans in by_file.items():
    spans.sort()
    for (a1, b1, f1), (a2, b2, f2) in zip(spans, spans[1:]):
        if a2 < b1:
            problems.append(f"overlap {f} {f1} {f2}")
print("problems:", problems)
for fid in E:
    for t in E[fid]:
        t.pop("line_found", None)
json.dump(E, open(os.path.join(HERE, "edits_repair.json"), "w"), indent=1, ensure_ascii=False)
