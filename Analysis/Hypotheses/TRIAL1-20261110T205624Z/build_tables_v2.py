"""Builds trial_outcomes_v2.tsv, excluded_v2.tsv and CHANGES_v2.tsv from the closed v1 tables.
Every new value was read on 2026-11-10 from a source saved under sources_v2/ (registry JSON, SAP/protocol
text, PMC full text, PubMed XML, ChEMBL JSON) or from the CTIS page saved in sources_v2/web/.
Row keys are 1-based data-row numbers of trial_outcomes.tsv (v1)."""
import csv, os
D = os.path.dirname(os.path.abspath(__file__))
NV = "not verified"
rows = list(csv.DictReader(open(os.path.join(D, "trial_outcomes.tsv")), delimiter="\t"))
cols = list(rows[0].keys())

S = dict(
    tirz_sap="NCT04166773 SAP v2 (sources_v2/docs) sections 3.4, 3.6; registry results NCT04166773",
    tirz_reg="registry results NCT04166773; PubMed 38856224",
    peg="registry results NCT04929483; SAP NCT04929483; PMC10718287 full text (PMID 37356033)",
    efx_abs="PubMed abstract PMID 37802088 (efruxifermin replicates FGF21 agonism of FGFR1c, 2c, 3c); ChEMBL CHEMBL4594389 has no mechanism",
    efx_sap="NCT05039450 SAP (no multiplicity adjustment)",
    resolveit="PubMed search NCT02704403 (only PMID 38061448, NIS2+ biomarker paper); NCT02704403 SAP revision notes",
    golden="registry results NCT01694849 (ALT LS mean difference, mixed model)",
    trop="registry results NCT02855164; SAP NCT02855164; PMC9941046 full text (PMID 36797481)",
    ervo="registry results NCT04321031; SAP NCT04321031 ('no formal hypothesis tests planned'); PubMed abstract PMID 40753985",
    deni="registry results NCT04906421; SAP NCT04906421 section 7.3",
    dean="PMC12135075 full text (PMID 40467095), fixed-sequence hierarchy",
    msdc="PubMed abstract PMID 31697972; registry results NCT02784444",
    bela="registry results NCT02462967",
    bi="PubMed 37932258; Nat Commun full text PMC10628239",
    chembl="ChEMBL MCP and REST name search 2026-11-10 (sources_v2/chembl)",
    fortuna="registry results NCT05809934; CTIS 2023-509704-14-00 (sources_v2/web); PubMed search NCT05809934 and AZD2693",
    aramchol="registry results NCT02279524; Prot_SAP NCT02279524; PMC12165723 full text (PMID 34621052)",
)

# (row, field, new value, source key)
CH = [
 # tirzepatide
 (8, "result_vs_control", "not_improved", "tirz_sap"),
 (8, "ci_or_p", "prespecified (evaluable participants, no imputation; no multiplicity adjustment): registry logistic 15 mg OR 2.47 (95% CI 0.95 to 6.40) p=0.063, 10 mg p=0.074, 5 mg OR 3.01 (1.15 to 7.90) p=0.025; paper (placebo-pattern imputation) 15 mg difference 95% CI 1 to 42", "tirz_sap"),
 (8, "notes", "key secondary; SAP v2 3.4 sets the evaluable-participant logistic analysis as primary and the imputed analysis as supplementary; coded on the highest dose (15 mg), not significant", "tirz_sap"),
 (9, "pmid", "38856224", "tirz_reg"),
 (9, "doi", "10.1056/NEJMoa2401943", "tirz_reg"),
 (9, "notes", "secondary endpoint; value from registry results, not stated in the abstract; SAP 3.6 no multiplicity adjustment", "tirz_sap"),
 # pegozafermin
 (14, "ci_or_p", "30 mg: 95% CI 9 to 33 (abstract); registry CMH with multiple imputation: 30 mg difference 20.89 (9.05 to 32.72) p=0.0009, 44 mg Q2W p=0.0005, 15 mg p<0.0001; SAP: each dose at two-sided 0.05, no multiplicity adjustment", "peg"),
 (14, "notes", "co-primary; full text reports P values for the fibrosis endpoint only and calls other estimates unadjusted", "peg"),
 (13, "notes", "ENLIVEN; registry 30 mg p=0.0085, 44 mg Q2W p=0.0077; ChEMBL has no mechanism record; no receptor named in abstract, full text or registry, so target_genes stays not verified", "peg"),
 # efruxifermin target
 (15, "target_genes", "FGFR1;FGFR2;FGFR3", "efx_abs"),
 (16, "target_genes", "FGFR1;FGFR2;FGFR3", "efx_abs"),
 (17, "target_genes", "FGFR1;FGFR2;FGFR3", "efx_abs"),
 (17, "notes", "SYMMETRY cohort D, add-on to GLP-1RA; secondary endpoint (primary was safety); SAP: no multiplicity adjustment; target from abstract PMID 37802088 (FGFR1c/2c/3c), no ChEMBL mechanism", "efx_sap"),
 # RESOLVE-IT
 (28, "notes", "RESOLVE-IT secondary endpoint; SAP revision: key secondaries gated at two-sided 0.01 after the primary, which failed; no primary results publication found in PubMed", "resolveit"),
 # GOLDEN-505 ALT
 (29, "result_vs_control", "not_improved", "golden"),
 (29, "effect_size", "120 mg LS mean difference in ALT change vs placebo -9.45 U/L (80 mg -6.13)", "golden"),
 (29, "ci_or_p", "120 mg 95% CI -19.4 to 0.49, p=0.062; 80 mg -15.97 to 3.71, p=0.221 (registry)", "golden"),
 (29, "notes", "GOLDEN-505; abstract claim of liver enzyme reduction is driven by ALP and GGT (registry p<0.001); ALT itself not significant", "golden"),
 # tropifexor
 (36, "result_vs_control", "improved_unadjusted", "trop"),
 (36, "ci_or_p", "registry: 200 ug vs placebo pooled parts A+B+C p=0.003 (part C p=0.030); unadjusted, full text calls p-values descriptive", "trop"),
 (37, "result_vs_control", "improved_unadjusted", "trop"),
 (37, "ci_or_p", "registry: 200 ug vs placebo part C p<0.001 (pooled p<0.001); unadjusted, full text calls p-values descriptive", "trop"),
 # ervogastat
 (40, "effect_size", "300 mg BID 4/31 vs placebo 12/34; risk difference -0.22 (favours placebo)", "ervo"),
 (40, "ci_or_p", "logistic 90% CI -0.30 to -0.05; Bayesian dose-response RD -0.10 (90% CrI -0.26 to 0.05); no formal hypothesis test planned", "ervo"),
 (40, "notes", "MIRNA; fibrosis component reported separately in registry; not coded worsened because no test was planned; DGAT2 target from abstract PMID 40753985; ChEMBL CHEMBL4760665 has no mechanism", "ervo"),
 (41, "result_vs_control", "improved_unadjusted", "ervo"),
 (41, "effect_size", "300 mg BID 13/31 vs placebo 3/34; risk difference 0.34", "ervo"),
 (41, "ci_or_p", "logistic 90% CI 0.10 to 0.61; no formal hypothesis test planned (SAP)", "ervo"),
 (41, "notes", "MIRNA secondary; 90% CI, not a significance test", "ervo"),
 # chembl
 (42, "chembl_id", "no ChEMBL record found (name search 2026-11-10)", "chembl"),
 (70, "chembl_id", "no ChEMBL record found (name search 2026-11-10)", "chembl"),
 (71, "chembl_id", "no ChEMBL record found (name search 2026-11-10)", "chembl"),
 (72, "chembl_id", "no ChEMBL record found (name search 2026-11-10)", "chembl"),
 # aramchol
 (44, "notes", "secondary; hierarchy stopped at primary (p=0.066); registry p=0.051, unadjusted", "aramchol"),
 (45, "notes", "secondary; hierarchy stopped at primary; registry p=0.21", "aramchol"),
 # denifanstat
 (48, "effect_size", "fibrosis improvement >=1 stage without worsening MASH, mITT: 33/81 vs 8/45", "deni"),
 (48, "ci_or_p", "two-sided CMH p=0.0102 (registry); SAP: only the co-primaries controlled for multiplicity", "deni"),
 (48, "notes", "FASCINATE-2 secondary", "deni"),
 # belapectin
 (57, "effect_size", ">=1-stage change in Brunt-Kleiner fibrosis: 8 mg/kg 2/54 vs placebo 1/54 (2 mg/kg 3/54)", "bela"),
 (57, "ci_or_p", "8 mg/kg p=0.602; 2 mg/kg p=0.362 (registry)", "bela"),
 # DEAN
 (59, "result_vs_control", "improved", "dean"),
 (59, "notes", "DEAN; fixed-sequence hierarchy primary -> resolution -> fibrosis at 0.05, all passed (full text)", "dean"),
 (60, "result_vs_control", "improved", "dean"),
 (60, "notes", "DEAN; third step of fixed-sequence hierarchy, passed (full text)", "dean"),
 # MSDC-0602K
 (63, "target_genes", "MPC1;MPC2", "msdc"),
 (64, "target_genes", "MPC1;MPC2", "msdc"),
 (65, "target_genes", "MPC1;MPC2", "msdc"),
 (63, "effect_size", "250 mg 27/86 vs placebo 15/74; OR 1.82", "msdc"),
 (63, "ci_or_p", "95% CI 0.86 to 3.82; p=0.116 (registry, no multiplicity adjustment)", "msdc"),
 (63, "notes", "EMMINENCE secondary; target from abstract ('designed to preferentially target the mitochondrial pyruvate carrier'); v1 note that the abstract names no target was wrong", "msdc"),
 (64, "effect_size", "250 mg 25/86 vs placebo 16/74; OR 1.51", "msdc"),
 (64, "ci_or_p", "95% CI 0.72 to 3.16; p=0.270 (registry, no multiplicity adjustment)", "msdc"),
 # BI 1467335
 (69, "pmid", "37932258", "bi"),
 (69, "doi", "10.1038/s41467-023-42398-w", "bi"),
 (69, "notes", "secondary endpoint; full text: 10 mg placebo-corrected relative ALT change -15.1% (90% CI -27.8 to -2.3); target AOC3 from registry description", "bi"),
 # FORTUNA
 (70, "notes", "FORTUNA primary; no publication in PubMed; registry other pre-specified lab outcome (baseline to 64 weeks): ALT >3x ULN 25 mg 18/73, 50 mg 33/69, placebo 10/74; AST >3x ULN 11/18/11; CTIS urgent safety measure (event 13 Sept 2025): transaminases persistently elevated in some participants at last safety visit, 14 identified for follow-up", "fortuna"),
]

FULL = {6, 13, 14, 20, 21, 22, 23, 24, 25, 26, 36, 37, 38, 39, 43, 44, 45, 46, 59, 60, 69}
REGR = {8, 9, 10, 11, 15, 16, 17, 18, 19, 27, 28, 29, 31, 35, 40, 41, 47, 48, 57, 63, 64, 65, 70, 71}
# full texts read: semaglutide phase 1, pegozafermin, PIVENS, tropifexor, firsocostat phase 2, aramchol, DEAN, BI 1467335
FULL -= {20, 21, 22, 38}   # pioglitazone 2016 and ATLAS full texts were not read
changes = []
for r, f, new, src in CH:
    old = rows[r - 1][f]
    if old != new:
        changes.append((r, rows[r - 1]["drug"], rows[r - 1]["endpoint"], f, old, new, S[src]))
        rows[r - 1][f] = new
for i, x in enumerate(rows, 1):
    x["verification_v2"] = "full text" if i in FULL else ("registry results" if i in REGR else "abstract")

with open(os.path.join(D, "trial_outcomes_v2.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, cols + ["verification_v2"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
with open(os.path.join(D, "CHANGES_v2.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["row_key", "drug", "endpoint", "field", "old", "new", "source_opened"]); w.writerows(changes)

# ---- excluded_v2
ex = list(csv.DictReader(open(os.path.join(D, "excluded.tsv")), delimiter="\t"))
excols = list(ex[0].keys())
upd = {
 "HSD17B13 knockdown": ("phase 1 to 2 only; no placebo-controlled readout on the four endpoints",
   "Rapirosiran phase 1 NCT04565717 (PMID 40581300): primary AEs, mRNA knockdown only; ARO-HSD NCT04202354: placebo only in healthy volunteers, NASH cohorts open label, results PK/safety; GSK4532990 HORIZON NCT05583344 active not recruiting, no results; SKYLINE NCT06104319 open label; AZD7503 NCT05864391 terminated (n=40), no results, NCT06093542 healthy volunteers, NCT05560607 open label; ALN-HSD NASHGEN-2 NCT05519475 recruiting; INI-822 NCT05945537 phase 1 recruiting"),
 "CIDEB knockdown": ("no randomized readout verified",
   "NCT06836609 phase 1/2 active, no posted results; NEJM correspondence PMID 43001761 has no abstract and no PMC full text, not opened; NCT07631637 recruiting"),
 "MTARC1 knockdown": ("no interventional drug trial found", "ClinicalTrials.gov searches (2026-11-10) for MTARC1 and drug names returned no drug trial"),
 "INHBE knockdown": ("no MASLD readout", "NCT06700538 phase 1/2 in obesity; NCT07859813 MASLD not yet recruiting; both opened 2026-11-10"),
 "GPAM knockdown": ("no interventional drug trial found", "ClinicalTrials.gov searches (2026-11-10) for GPAM and drug names returned no drug trial"),
 "MOGAT2 inhibitor": ("phase 1 only; no placebo-controlled MASLD readout",
   "BMS-963272: NCT04116632, NCT04124003, NCT02327273 phase 1; NCT04766476 phase 1b NAFLD randomized, terminated at n=9 ('Business objectives have changed'), no results; PMID 36323235 reports phase 1 obesity data only"),
 "RORA": ("phase 1 healthy volunteers only", "TB-840 (RORalpha agonist) phase 1 NCT05045534, PMID 41010352; no MASLD trial"),
 "CB1 (CNR1)": ("randomized NASH trials terminated without results",
   "rimonabant phase 3 NASH NCT00576667 (n=165) and NCT00577148 (n=89) terminated ('Company decision' after health authority action), no posted results, no PubMed hit by NCT; PMID 22828945 is an obesity RCT with MRS liver fat, not a NAFLD-selected population; RYI-018 NCT03261739 phase 1, no results"),
}
for x in ex:
    if x["mechanism"] in upd:
        x["status"], x[excols[4]] = upd[x["mechanism"]]
with open(os.path.join(D, "excluded_v2.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, excols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(ex)

nv = sum(v == NV for x in rows for k, v in x.items())
print("rows", len(rows), "changes", len(changes), "not verified fields", nv)
from collections import Counter
print(Counter(x["result_vs_control"] for x in rows))
print([(i, k) for i, x in enumerate(rows, 1) for k, v in x.items() if v == NV])
