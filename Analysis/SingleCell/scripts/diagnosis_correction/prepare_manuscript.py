"""Prepare reviewable manuscript corrections without changing the named release."""
import hashlib
import json
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[4]
task = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
d = json.loads((task / "corrected_results.json").read_text())
out = task / "manuscript"
out.mkdir(exist_ok=False)
changes = []
hashes = []


def read(name):
    p = root / "docs/manuscript" / name
    hashes.append(dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    return p.read_text()


def paragraph(text, prefix, replacement, file):
    parts = text.split("\n\n")
    indices = [i for i, x in enumerate(parts) if x.startswith(prefix)]
    assert len(indices) == 1, (file, prefix)
    i = indices[0]
    changes.append(dict(file=file, before=parts[i], after=replacement))
    parts[i] = replacement
    return "\n\n".join(parts)


def save(name, text):
    p = out / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


def focal(rows, module):
    return next(x for x in rows if x["module"] == module)


def effect(r):
    return f"beta = {r['gls_beta']:.3f} (nominal 95% CI {r['gls_ci_low']:.3f} to {r['gls_ci_high']:.3f}; permutation p = {r['perm_p']:.5g}, BH q = {r['perm_q']:.3f})"


p = d["programs"]; e = focal(p["focal"], 8); b = focal(p["focal"], 20)
q = d["qc"]; a = d["ambient"]; t = d["tcell20"]
names = {"hepatocytes": "hepatocyte", "cholangiocytes": "cholangiocyte", "fibroblasts": "fibroblast", "macrophages": "macrophage", "tcells": "T-cell"}
ns = ", ".join(f"{r['n_fit']} {names[r['cell_type']]}" for r in d["primary_n"])
stages = ", ".join(f"{r['n']} {r['disease_stage_coarse'].replace('_', ' ')}" for r in d["descriptive_diagnoses"])
source_note = (
    "Deposited diagnosis, measured histology and etiology were kept separate. GSE202379 P98 was excluded from disease inference because its Healthy-control label conflicts with SAF S1A1F1. "
    "P30 (S0A1F0) was retained as a qualified nonsteatotic reference, with an exclusion sensitivity; it is not confirmed histologically normal tissue. "
    "Official RNA/GSM metadata place local GSE244832 D11/JB_321 in MASL and D14/JB_337 in MASH. This correction does not verify RNA–ATAC pair identities. "
    "All 19 GSE192742 donors remain outside disease inference. Alcohol/PBC cirrhosis is not classified as MASLD, and GSE189600 replicate independence remains unresolved."
)
primary_text = (
    "Two focal programs, Stromal ECM (IGFBP7) and Ductular injury (BICC1), were discovered in hepatocyte profiles but are among the 37 whose content contradicts that label. "
    "We tested all 117 fixed programs for cross-sectional association with Healthy/Steatosis/Steatohepatitis categories after correcting source diagnoses. "
    f"The ambient-correctable universe contains at most 59 source-qualified donors across four datasets; the two multi-stage datasets, GSE202379 and GSE244832, supply the fitted slopes ({ns} donors). "
    "We fitted score ~ ordinal stage + dataset by generalized least squares (REML; nlme::gls with varIdent; Pinheiro et al.), with a separate residual variance per dataset and 19,999 within-dataset stage permutations (seed 20260923). "
    f"No program passes BH correction over 117 programs plus the unchanged stability rule (smallest q = {p['smallest_q']:.3f}). Ductular injury has {effect(b)}; Stromal ECM has {effect(e)} (Figure 4C). "
    f"The two cohort slopes have opposite signs for {p['cohort_opposite_signs']} programs; {p['cohort_Q_BH05']} heterogeneity tests pass BH over 117 (smallest q = {p['cohort_smallest_Q_q']:.3f}). "
    "The variance model was chosen after earlier focal-program slopes and equal-variance results were known. These corrected estimates remain conditional on two cohorts, their source diagnoses and the ordinal coding. Failure to pass correction does not establish absent effects."
)
ambient_parts = []
for module, label in [(8, "Stromal ECM"), (20, "Ductular injury")]:
    r = focal(d["ambient_paired_focal"], module); z = focal(a["focal"], module)
    ambient_parts.append(f"{label} changes from {r['raw_beta']:.4f} to {r['corrected_beta']:.4f} in the paired equal-variance analysis (delta {r['delta_beta']:.4f}, nominal p = {r['delta_pvalue']:.4g}, BH q = {r['delta_qvalue']:.3f}); its ambient-corrected GLS beta is {z['gls_beta']:.3f} (permutation BH q = {z['perm_q']:.3f}).")
ambient_text = (
    "We next examined ambient RNA using unchanged scores, original full-transcriptome exposures and fixed program weights in the correction-available donor universe. "
    + " ".join(ambient_parts) +
    f" The GLS sensitivity tests 106 non-T-cell programs within the 117-program family and selects {a['selected']} (smallest q = {a['smallest_q']:.3f}). "
    "The 11 T-cell programs are untestable on the corrected native scale. Attenuation classification is conditional on primary selection, so it applies to neither focal program. "
    "These results do not establish the cellular origin of either program's transcripts."
)
lineage_parts = []
for r in d["content_lineage"]:
    label = "Stromal ECM in fibroblasts" if "f05c535" in r["program_uid"] else "Ductular injury in cholangiocytes"
    lineage_parts.append(f"{label}: beta {r['beta']:.3f} (HC3-based meta-analysis 95% CI {r['low']:.3f} to {r['high']:.3f}; n = {r['n']} donors; two-program BH q = {r['q_two_programs']:.3f}).")
lineage_text = (
    "We compared steatohepatitis with steatosis within each focal program's content-associated lineage, using the same correction-available cells and donors and the corrected RNA diagnoses. "
    + " ".join(lineage_parts) +
    " These are targeted projected-expression tests in GSE202379 and GSE244832, with a 50-cell primary gate; they do not independently replicate the hepatocyte discovery or exclude within-lineage substate mixtures. "
    "Cohort-specific estimates, full-raw, 20-cell and annotation-confidence sensitivities are reported separately. Captured-lineage fractions are not whole-tissue proportions."
)
c = d["communication"]; supported = c["supported"]; assert len(supported) == 1
r = supported[0]
ccc_text = (
    f"A supplementary donor-level communication screen contains {c['native_groups']:,} native ligand-receptor/cell-type groups: {c['native_fitted']} are fitted and {c['native_groups']-c['native_fitted']:,} have insufficient stage counts. "
    "It retains the native within-cell-type-pair/term BH families followed by across-cell-pair Bonferroni correction. "
    f"The separate fixed 26-test program screen has {c['estimable']} estimable tests. One Ductular injury association passes model-based and HC3 BH: hepatocyte CDH1–endothelial PTPRM, beta = {r['beta']:.4f}, n = {r['n_donors']}, q = {r['qvalue']:.4g}, HC3 q = {r['hc3_qvalue']:.4g}. "
    f"Its sign agrees in {r['n_lodo_same_direction']} of {r['n_lodo_estimable']} cohort omissions (coefficient range {r['lodo_beta_min']:.4f} to {r['lodo_beta_max']:.4f}), so cohort stability is limited. "
    "These same-atlas RNA associations do not establish ligand secretion, receptor activation or signaling. Of four selected transcription-factor regulons, two remain untestable because of insufficient retained weight."
)
file = "draft/fig4.md"; text = read(file)
old = next(x for x in text.split("\n\n") if x.startswith("To place disease-associated"))
text = paragraph(text, "To place disease-associated", old[:old.index("Primary-source diagnoses")] + f"Source-qualified descriptive categories are {stages}. " + source_note, file)
text = paragraph(text, "Two focal programs, Stromal ECM", primary_text, file)
text = paragraph(text, "We next asked how much ambient", ambient_text, file)
text = paragraph(text, "We next compared steatohepatitis", lineage_text, file)
text = paragraph(text, "A supplementary cell-cell communication", ccc_text, file)
save(file, text)
file = "draft/01_intro.md"; text = read(file)
assert "single-cell data from 60 donors" in text
text = text.replace("single-cell data from 60 donors", "single-cell data from at most 59 source-qualified donors")
save(file, text)
file = "draft/supplementary_notes.md"; text = read(file)
text = text.replace("smallest q of any program 0.323", f"smallest q of any program {p['smallest_q']:.3f}")
text = paragraph(text, "Units in these analyses are biological donors.", "Units in the primary analyses are source-qualified donor identifiers. The 61-record equal-variance sensitivity includes two GSE189600 records whose biological independence is unresolved.", file)
text = paragraph(text, "The 62-donor equal-variance fit", "The corrected equal-variance sensitivity retains at most 61 source records, including two GSE189600 records without ambient correction and with unresolved independence. Its complete 117-program ordinary and HC3 estimates are reported separately. Historical 62-record p/q values and cell-gate comparisons are superseded; they do not supply current primary evidence.", file)
text = paragraph(text, "Under the primary model", f"Under the corrected primary model, averaging stored scores over QC-passing cells selects {q['selected']} of 117 programs (smallest permutation BH q = {q['smallest_q']:.3f}); stage betas correlate with the primary estimates at Spearman {q['beta_rho']:.3f}. Ductular injury has {effect(focal(q['focal'],20))}; Stromal ECM has {effect(focal(q['focal'],8))}. These stored cell scores include neighbor contributions from failing cells, so this is not rescoring after QC. Removing donors with fewer than 20 T cells leaves {t['n']} of 52 fitted T-cell donors, with {t['selected']} selected T-cell programs (smallest T-cell q after BH over the complete 117-program family = {t['smallest_q']:.3f}); the other 106 permutation p values remain fixed.", file)
text = paragraph(text, "The Stromal ECM paired change", ambient_text, file)
text = paragraph(text, "Both content-lineage results", lineage_text, file)
save(file, text)
file = "METHODS.md"; text = read(file)
text = paragraph(text, "Seven public human liver single-cell", "Seven public human liver single-cell or single-nucleus datasets contribute 275 libraries, 1,232,318 cells and 104 source-based donor identifiers. Two low-cell GSE189600 libraries are excluded from the analyzed atlas, leaving 273 libraries, 1,232,285 cells and 102 identifiers. GSE189600 replicate independence is unresolved; its records are excluded from primary disease, stage and composition inference and are retained only in the explicitly qualified 61-record equal-variance sensitivity.", file)
old = next(x for x in text.split("\n\n") if x.startswith("Program discovery and scoring used"))
text = paragraph(text, "Program discovery and scoring used", old[:old.index("Averaging the stored")] + "QC-passing-cell averages are reanalyzed under the corrected primary GLS model (Supplementary Note 8). Stored cell scores include neighbor contributions from failing cells, so this is not rescoring after QC.", file)
text = paragraph(text, "Several datasets deposit multiple libraries", "Libraries were joined to source-based donor identifiers before analysis. Eligible preparation libraries were filtered before their counts were summed within donors, so excluded sort fractions could not reenter the count matrix during collapse. Per-cell-type differential expression used edgeR normalization and limma-voom donor pseudobulk with dataset and disease covariates; BH correction was applied over all tested genes separately for each cell type and declared contrast. Designs required estimable groups and full rank. The hepatocyte binary comparison contains 73 donors (15 qualified references and 58 MASLD cases); source nonsteatotic P30 is excluded in a separate sensitivity. Whole-liver composition tests exclude sorted fractions, GSE192742, GSE136103 and GSE189600; 70 donors (13 references, 57 cases) enter the dataset-adjusted arcsine-square-root proportion model with limma empirical Bayes. BH families are the 16 lineages, frozen hepatocyte clusters and frozen metaclusters separately. These are captured-cell fractions, not measured tissue composition. Statistical significance in one lineage does not establish cell-type specificity.", file)
text = paragraph(text, "We corrected diagnosis labels against", f"The source-qualified descriptive categories among the 102 analyzed identifiers are {stages}. " + source_note + " The GSE174748 disease explants remain source-documented Cirrhosis and outside the primary three-category estimand.", file)
text = paragraph(text, "F-stage-specific single-cell analyses", "Measured fibrosis analyses use only deposited numeric SAF entries in GSE202379: 39 source donors. Free-text healthy-control entries for PHL1/PHL2 and five free-text end-stage entries do not establish measured F0 or F4. P98 has measured F1 but remains excluded from disease inference because of its diagnosis conflict. Twelve source donors have measured F0/F1, below the minimum of 15 for H7. Numeric-SAF sensitivities use their own eligible donor rosters. Predicted or disease-augmented F stages are excluded from primary cross-cohort inference.", file)
text = paragraph(text, "The primary ordinal models use", primary_text + " Selection retains the frozen stability median of at least 0.5 (best-match Jaccard of the top 50 genes when each dataset is omitted). Each slope cohort must have at least eight donors and at least two stage categories containing two or more donors. Ineligible programs retain p = 1 in the 117-program family.", file)
nc = pd.read_csv(task / "programs_diagnostics/NC1_null_family.tsv", sep="\t").iloc[0]
text = paragraph(text, "Before the outcome was read", f"Using the corrected donor labels, the first 2,000 complete null draws gave at least one 117-program BH call in {100*nc.fraction_any_call:.2f}% of draws. Known-truth simulation retained 2,000 draws for each of beta 0 and beta 0.30; model-based 95% coverage was 94.4% and 94.25%, respectively, with no failed fits. Detection experiments retain the previously fixed alternatives, real scores, 1,000 within-dataset draws and the permutation threshold at 0.05/117. They measure conditional sensitivity, not absence of effects. Seeds are 20260923. The two focal programs were selected in an earlier analysis and were not reselected after correction.", file)
text = paragraph(text, "In the 62-donor equal-variance sensitivity", "Primary-model sensitivities average stored scores over QC-passing cells or exclude donors with fewer than 20 T cells. Each preserves the 117-program BH family; the T-cell check changes only the 11 T-cell permutation p values. The corrected equal-variance sensitivity includes 61 source records, including two GSE189600 records of unresolved independence. Historical gate-specific significance claims are not carried forward. Projected scores are not native Hotspot units.", file)
old = next(x for x in text.split("\n\n") if x.startswith("The atlas used for Hotspot discovery"))
old = old.replace("at most 60 donors", "at most 59 source-qualified donors")
old = old.replace("the eligibility rule removes only the two source-cirrhotic donors", "eligibility also excludes P98 and follows source diagnosis, etiology and dataset restrictions")
text = paragraph(text, "The atlas used for Hotspot discovery", old, file)
text = paragraph(text, "We inferred cell-cell communication", ccc_text, file)
text = paragraph(text, "The screen has 22 estimable", "The corrected fixed screen has 20 estimable and six untestable tests in the 26-test family. The one supported association and its limited leave-one-cohort-out stability are reported with the complete family in Figure S4. Untestable associations are not negative evidence.", file)
save(file, text)
file = "05_figure_legends.md"; text = read(file)
prefix = "**(C)** Ordinal stage β for the 117 programs"
text = paragraph(text, prefix, f"**(C)** Twenty-four fixed display programs from the complete 117-program family. GLS ordinal stage beta with per-dataset residual variance, using source-qualified diagnoses in GSE202379 and GSE244832 (54 hepatocyte donors; other lineage n in source table). Gray bars use the within-dataset permutation critical |t| at 0.05/117; black bars are nominal 95% CIs. No program passes the complete-family permutation BH and stability rule (smallest q = {p['smallest_q']:.3f}); all 117 effects are shown in S4C. The fewer-than-20-T-cell sensitivity retains {t['n']} donors and selects no T-cell program (smallest T-cell q after BH over 117 = {t['smallest_q']:.3f}).", file)
text = paragraph(text, "**(A)** Hepatocyte native-score effect", "**(A)** Corrected equal-variance sensitivity for Stromal ECM and Ductular injury, with at most 61 source records. Two GSE189600 records have unresolved biological independence. Pooled, cohort and cohort-omission estimates retain their actual denominators and nominal HC3 95% intervals; these are not complete-family primary calls.", file)
text = paragraph(text, "**(B)** MASH minus steatosis projected score", "**(B)** Corrected MASH-minus-steatosis projected scores in fibroblasts (Stromal ECM) and cholangiocytes (Ductular injury), using the same correction-available cells and donors at the 50-cell gate. Cohort-specific HC3 95% CIs and full donor rosters accompany the source table. The two programs form one BH family per declared arm. Projected scores are not native Hotspot units.", file)
text = paragraph(text, "**(C)** Median dataset-centered donor score", "**(C)** Donor-score means after centering each program across the fitted donor universe, shown for Healthy, Steatosis and Steatohepatitis source categories for all 117 fixed programs. The heatmap is descriptive; the accompanying GLS stage beta adjusts for dataset with separate residual variances. Black bars, nominal 95% CI; gray bars, 117-program permutation interval. Source metadata and fitted denominators match Figure 4C.", file)
text = paragraph(text, "**(E)** Adjusted associations of the two focal", "**(E)** Corrected program–communication associations in the 26-test BH family, with 20 estimable and six untestable tests. One Ductular injury association passes both ordinary and HC3 BH: CDH1–PTPRM, 57 donors, beta 0.4984, q 9.940e-6, HC3 q 5.952e-10. Its sign agrees in three of four cohort omissions. Bars are nominal model-based 95% CIs; untestable rows are listed in the complete source table. RNA covariation does not establish signaling.", file)
save(file, text)
pd.DataFrame(changes).to_csv(out / "paragraph_changes.tsv", sep="\t", index=False)
pd.DataFrame(hashes).to_csv(out / "source_document_hashes.tsv", sep="\t", index=False)
(out / "REVIEW_STATE.txt").write_text("Corrected text candidate only. Current named release and working manuscript files were not replaced. Final adoption requires explicit user approval.\n")
print("Prepared corrected manuscript paragraphs and matching figure legends for review.")
