#!/usr/bin/env python3
"""H26: therapeutic direction calls for genetically anchored candidates (provisional on the current
COLOC release). Implements PRESPEC_H26.md. Run with `python -I` (reads downloaded, untrusted tables).

Usage: python -I h26_direction_calls.py <result_dir> <new_run_dir>
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = REPO / "GWAS/finemapping"
CATALOG = REPO / "RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-20260925T175215Z/passport_gene_index.parquet"
COLOC = FM / "results/susie_coloc_r2_20260924T053244Z"
DIRECTION = REPO / ("Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/"
                    "coloc_direction/full-20260930T120450Z")
BROADAWAY = REPO / "data/broadaway_eqtl"
RELEASE = "candidate_r2_20260924T053244Z"
PROVISIONAL = "provisional on the current COLOC release (susie_coloc_r2_20260924T053244Z)"
TRAIT_PRIORITY = ["NAFLD", "NASH", "PDFF", "ALT", "AST", "GGT"]
CODING = {"PNPLA3": "I148M", "TM6SF2": "E167K", "MTARC1": "A165T", "GPAM": "I43V",
          "HSD17B13": "rs72613567 splice", "GCKR": "P446L", "SERPINA1": "Z allele", "APOE": "isoforms"}
CALIBRATION = {  # gene: expected sign of pLoF burden on ALT / liver fat (None = no expectation)
    "PNPLA3": None, "HSD17B13": -1, "MTARC1": -1, "GPAM": -1, "APOB": +1, "CIDEB": -1}
CHEN_COHORTS = ["gb_ukbb_pdff", "gb_ukbb_masld", "gb_aou_masld", "gb_regeneron_masld", "gb_sema4_masld"]
# Recalled public knowledge; every row is unverified until checked against ChEMBL / ClinicalTrials.gov.
# direction: inhibit = antagonist / inhibitor / silencing; activate = agonist / analogue.
DRUG_TARGETS = [
    ("THRB", "resmetirom; VK2809", "activate", "approved / phase 2-3"),
    ("GLP1R", "semaglutide; tirzepatide; survodutide", "activate", "approved / phase 3"),
    ("GIPR", "tirzepatide", "activate", "phase 3"),
    ("GCGR", "survodutide; pemvidutide", "activate", "phase 3 / phase 2"),
    ("KLB", "efruxifermin; pegozafermin (FGF21 analogues, FGFR1c/KLB)", "activate", "phase 3"),
    ("FGFR1", "efruxifermin; pegozafermin", "activate", "phase 3"),
    ("FGFR4", "aldafermin (FGF19 analogue)", "activate", "phase 2b did not meet fibrosis endpoint"),
    ("PPARA", "lanifibranor; elafibranor", "activate", "phase 3 / failed phase 3 (elafibranor)"),
    ("PPARD", "lanifibranor; elafibranor", "activate", "phase 3 / failed phase 3 (elafibranor)"),
    ("PPARG", "lanifibranor; pioglitazone", "activate", "phase 3"),
    ("NR1H4", "obeticholic acid; cilofexor; tropifexor", "activate", "not approved for NASH / phase 2"),
    ("ACACA", "firsocostat", "inhibit", "phase 2"),
    ("ACACB", "firsocostat", "inhibit", "phase 2"),
    ("DGAT2", "ervogastat; ION224", "inhibit", "phase 2"),
    ("SCD", "aramchol", "inhibit", "phase 3 stopped"),
    ("FASN", "denifanstat", "inhibit", "phase 3"),
    ("KHK", "PF-06835919", "inhibit", "phase 2"),
    ("MAP3K5", "selonsertib", "inhibit", "failed phase 3"),
    ("CCR2", "cenicriviroc", "inhibit", "failed phase 3"),
    ("CCR5", "cenicriviroc", "inhibit", "failed phase 3"),
    ("LGALS3", "belapectin", "inhibit", "failed phase 2b/3"),
    ("LOXL2", "simtuzumab", "inhibit", "failed phase 2b"),
    ("AOC3", "BI 1467335", "inhibit", "discontinued phase 2"),
    ("MPC1", "MSDC-0602K", "inhibit", "phase 2b did not meet primary endpoint"),
    ("MPC2", "MSDC-0602K", "inhibit", "phase 2b did not meet primary endpoint"),
    ("MOGAT2", "BMS-963272", "inhibit", "phase 1-2"),
    ("HSD17B13", "rapirosiran (ALN-HSD); GSK4532990 (ARO-HSD); INI-822", "inhibit", "phase 1-2"),
    ("PNPLA3", "AZD2693; ARO-PNPLA3 (JNJ-75220795)", "inhibit", "phase 1-2"),
    ("MTARC1", "siRNA programmes", "inhibit", "phase 1"),
    ("CIDEB", "ALN-CIDE", "inhibit", "phase 1"),
]


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = ok.sum()
    if n == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def sign(x):
    return 0 if not np.isfinite(x) or x == 0 else (1 if x > 0 else -1)


def call_from_expression_sign(s):   # trait-raising haplotype increases expression -> inhibit
    return {1: "inhibit", -1: "activate"}.get(s, "")


def call_from_burden_sign(s):       # burden-associated increase in the trait -> activate (do not inhibit)
    return {1: "activate", -1: "inhibit"}.get(s, "")


# ----------------------------------------------------------------------------- raw-file readers
def gwas_rows(study, positions, registry):
    """Rows of a raw GWAS file at the given (chr, pos) set."""
    path = FM / registry.loc[study, "sumstats_path"]
    want = {(str(c), int(p)) for c, p in positions}
    wanted_pos = {p for _, p in want}
    keep = []
    for chunk in pd.read_csv(path, sep="\t", dtype={"chromosome": str}, chunksize=2_000_000,
                             usecols=["chromosome", "position", "allele1", "allele2", "beta", "pval"]):
        chunk = chunk[chunk["position"].isin(wanted_pos)]
        if len(chunk):
            chunk = chunk[[(c, int(p)) in want for c, p in zip(chunk["chromosome"], chunk["position"])]]
            keep.append(chunk)
    cols = ["chromosome", "position", "allele1", "allele2", "beta", "pval"]
    return pd.concat(keep) if keep else pd.DataFrame(columns=cols)


def broadaway_rows(wanted):
    """wanted: set of (chr, pos, ensg). Returns the matching Broadaway marginal rows."""
    out = []
    by_chr = {}
    for c, p, g in wanted:
        by_chr.setdefault(str(c), set()).add((int(p), g))
    for c, items in by_chr.items():
        path = BROADAWAY / f"chr{c}_marginal_summary_results.tsv"
        if not path.exists():
            continue
        pos = {p for p, _ in items}
        for chunk in pd.read_csv(path, sep="\t", chunksize=2_000_000,
                                 usecols=["CHR", "POS", "NEA", "EA", "Beta", "PVAL", "GeneSymbol", "ENSG"]):
            chunk = chunk[chunk["POS"].isin(pos)]
            if len(chunk):
                out.append(chunk[[(int(p), g) in items for p, g in zip(chunk["POS"], chunk["ENSG"])]])
    return pd.concat(out) if out else pd.DataFrame()


# ----------------------------------------------------------------------------- burden readers
def read_backman(path):
    d = pd.read_csv(path, sep="\t", usecols=["effect_allele", "beta", "standard_error", "effect_allele_frequency",
                                             "p_value", "Name"])
    ext = d["Name"].str.extract(r"^(?P<symbol>.+)\((?P<ensembl>ENSG\d+)\)\.GENE\.(?P<mask>M\d\.\w+)$")
    d = pd.concat([d, ext], axis=1)
    if d["ensembl"].isna().any() or (d["mask"] != d["effect_allele"]).any():
        raise SystemExit(f"unexpected Name format in {path}")
    return d


def read_chen(zdir):
    per = []
    for c in CHEN_COHORTS:
        d = pd.read_csv(zdir / f"{c}.tsv", sep="\t")
        d = d[d["test"] == "ADD"].copy()
        d["z"] = d["beta"] / d["standard_error"]
        d = d[np.isfinite(d["z"])]
        if d["gene"].duplicated().any():
            raise SystemExit(f"duplicated ADD gene rows in {c}")
        d["cohort"] = c
        per.append(d[["gene", "cohort", "beta", "standard_error", "z", "n", "minor_allele_count"]])
    per = pd.concat(per)
    per["w"] = np.sqrt(per["n"])
    per["wz"], per["w2"] = per["w"] * per["z"], per["w"] ** 2
    g = per.groupby("gene").agg(wz=("wz", "sum"), w2=("w2", "sum"), n_cohorts=("cohort", "size"),
                                n_total=("n", "sum"), mac_total=("minor_allele_count", "sum"),
                                cohort_signs=("z", lambda z: "".join("+" if v > 0 else "-" for v in z)))
    g["chen_z"] = g["wz"] / np.sqrt(g["w2"])
    g["chen_p"] = 2 * stats.norm.sf(np.abs(g["chen_z"]))
    authors = pd.read_csv(zdir / "gb_meta_analysis_sample_size_doe.TBL", sep="\t").set_index("MarkerName")
    g["chen_authors_meta_z"] = authors["Zscore"].reindex(g.index)
    g["chen_authors_meta_p"] = authors["P-value"].reindex(g.index)
    return g.drop(columns=["wz", "w2"]).reset_index().rename(columns={"gene": "symbol"}), per


def main():
    root, out = Path(sys.argv[1]), Path(sys.argv[2])   # downloads live in root; outputs go to a new run directory
    zdir, bdir = root / "zenodo_14804238", root / "gwascatalog_backman2021"
    summary = {"label": PROVISIONAL, "unit": "gene"}

    cat = pd.read_parquet(CATALOG)[["ensembl_id", "symbol", "primary_evidence_class"]]
    anchored = cat[cat["primary_evidence_class"] == "genetically_anchored"].copy()
    concordant = cat[cat["primary_evidence_class"] == "concordant"].copy()
    state = cat[cat["primary_evidence_class"] == "established_state_associated"].copy()
    assert (len(anchored), len(concordant), len(state)) == (383, 23, 1054), (len(anchored), len(concordant), len(state))
    genes = pd.concat([anchored, concordant]).rename(columns={"ensembl_id": "ensembl"})

    tier = pd.read_csv(COLOC / "gene_level_coloc_tier12.csv")
    genes = genes.merge(tier[["ensembl", "coloc_best_susie_pp4", "driving_trait", "susie_driving_gwas", "susie_state_t12"]],
                        on="ensembl", how="left")
    summary["anchored_genes_susie_state"] = genes[genes["primary_evidence_class"] == "genetically_anchored"][
        "susie_state_t12"].value_counts(dropna=False).to_dict()

    registry = pd.read_csv(FM / "config/gwas_registry.tsv", sep="\t").set_index("study_name")
    audit = pd.read_csv(DIRECTION / "gwas_orientation_audit.tsv", sep="\t").set_index("study")

    # ------------------------------------------------------------------ harmonisation checks 1-3
    checks = []
    known = [("PNPLA3 rs738409", "22", 44324727, "G", "C", +1, ["UKBB_ALT", "MVP_ALT_EUR", "MVP_NAFLD_EUR"]),
             ("HSD17B13 rs72613567", "4", 88231392, "TA", "T", -1, ["UKBB_ALT", "MVP_ALT_EUR"])]
    studies_needed = {s for k in known for s in k[6]}
    pairs = pd.read_csv(DIRECTION / "direction_signal_pairs.tsv", sep="\t", low_memory=False)
    pairs = pairs[(pairs["release"] == RELEASE)]
    same = pairs[(pairs["same_lead"].astype(str) == "True") & pairs["direction_sign"].isin([1, -1])].copy()
    same[["c", "p"]] = same["hit1"].str.split(":", expand=True)
    pos_by_study = {s: set() for s in studies_needed | set(same["study"])}
    for k in known:
        for s in k[6]:
            pos_by_study[s].add((k[1], k[2]))
    for s, c, p in zip(same["study"], same["c"], same["p"]):
        pos_by_study[s].add((c, int(p)))
    raw = {s: gwas_rows(s, pos, registry) for s, pos in pos_by_study.items()}

    def audited_column(study):
        return audit.loc[study, "effect_allele_column"] if study in audit.index else "unresolved"

    for label, c, p, raising, other, expected, studies in known:
        for s in studies:
            r = raw[s]
            r = r[(r["chromosome"] == c) & (r["position"] == p)]
            r = r[[{a, b} == {raising, other} for a, b in zip(r["allele1"], r["allele2"])]]
            col = audited_column(s)
            if r.empty or col not in ("allele1", "allele2"):
                checks.append(dict(check=label, study=s, status="not available (variant absent or study unresolved)"))
                continue
            r = r.sort_values("pval").iloc[0]
            beta_raising = r["beta"] if r[col] == raising else -r["beta"]
            checks.append(dict(check=label, study=s, allele1=r["allele1"], allele2=r["allele2"], beta=r["beta"],
                               p=r["pval"], effect_allele_column=col, tested_allele=raising,
                               beta_tested_allele=beta_raising, expected_sign=expected,
                               status="pass" if sign(beta_raising) == expected else "FAIL"))

    wanted_eqtl = {(c, int(p), g) for c, p, g in zip(same["c"], same["p"], same["ensembl"])}
    sort1_id = cat.loc[cat["symbol"] == "SORT1", "ensembl_id"].iloc[0]
    wanted_eqtl.add(("1", 109817590, sort1_id))
    eq = broadaway_rows(wanted_eqtl)
    r = eq[(eq["ENSG"] == sort1_id) & (eq["POS"] == 109817590)]
    if r.empty:
        checks.append(dict(check="SORT1 rs12740374 eQTL", study="Broadaway", status="not available"))
    else:
        r = r.iloc[0]
        b = r["Beta"] if r["EA"] == "T" else -r["Beta"]
        checks.append(dict(check="SORT1 rs12740374 eQTL", study="Broadaway", allele1=r["EA"], allele2=r["NEA"],
                           beta=r["Beta"], p=r["PVAL"], tested_allele="T", beta_tested_allele=b, expected_sign=1,
                           status="pass" if (sign(b) == 1 and {r["EA"], r["NEA"]} == {"G", "T"}) else "FAIL"))
    checks = pd.DataFrame(checks)
    checks.to_csv(out / "harmonisation_checks.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ check 4: same-lead recompute
    rec = []
    for _, row in same.iterrows():
        g = raw[row["study"]]
        g = g[(g["chromosome"] == row["c"]) & (g["position"] == int(row["p"]))]
        e = eq[(eq["ENSG"] == row["ensembl"]) & (eq["POS"] == int(row["p"]))]
        col = audited_column(row["study"])
        best = None
        for _, gr in g.sort_values("pval").iterrows():
            m = e[[{a, b} == {gr["allele1"], gr["allele2"]} for a, b in zip(e["EA"], e["NEA"])]]
            if len(m):
                best = (gr, m.iloc[0])
                break
        if best is None or col not in ("allele1", "allele2"):
            rec.append(dict(gene=row["gene"], study=row["study"], hit=row["hit1"], table_sign=row["direction_sign"],
                            recomputed_sign=np.nan, status="not recomputable"))
            continue
        gr, er = best
        effect_allele = gr[col]
        eqtl_beta = er["Beta"] if er["EA"] == effect_allele else -er["Beta"]
        s = sign(gr["beta"]) * sign(eqtl_beta)
        rec.append(dict(gene=row["gene"], study=row["study"], hit=row["hit1"], gwas_allele1=gr["allele1"],
                        gwas_allele2=gr["allele2"], gwas_beta=gr["beta"], effect_allele=effect_allele,
                        eqtl_ea=er["EA"], eqtl_beta=er["Beta"], table_sign=row["direction_sign"], recomputed_sign=s,
                        status="match" if s == row["direction_sign"] else "MISMATCH"))
    rec = pd.DataFrame(rec)
    rec.to_csv(out / "same_lead_recheck.tsv", sep="\t", index=False)
    n_rec = int((rec["status"] != "not recomputable").sum())
    n_mis = int((rec["status"] == "MISMATCH").sum())
    summary["harmonisation"] = {
        "checks_1_to_3": checks[["check", "study", "status"]].to_dict("records"),
        "same_lead_pairs_with_call": len(rec), "recomputed": n_rec, "mismatch": n_mis,
        "line_A_withheld": bool((checks["status"] == "FAIL").any() or (n_rec and n_mis / n_rec > 0.02))}

    # ------------------------------------------------------------------ line A
    dt = pd.read_csv(DIRECTION / "direction_by_trait.tsv", sep="\t", low_memory=False)
    dt = dt[dt["release"] == RELEASE]
    rows = []
    for ens, d in dt.groupby("ensembl"):
        dd = d[d["direction_sign"].isin([1, -1])]
        signs = sorted(set(dd["direction_sign"].astype(int)))
        gtex_opposes = bool((dd["check_gtex_liver_same_sign_hit2"] == 0).any())
        if dd.empty:
            status, s = "untestable", 0
        elif len(signs) > 1:
            status, s = "ambiguous (traits disagree)", 0
        elif gtex_opposes:
            status, s = "ambiguous (GTEx liver eQTL opposite sign)", 0
        else:
            status, s = "call", signs[0]
        traits = [t for t in TRAIT_PRIORITY if t in set(dd["trait"])]
        rows.append(dict(ensembl=ens, lineA_status=status, lineA_sign=s, lineA_call=call_from_expression_sign(s),
                         lineA_defining_trait=traits[0] if traits and status == "call" else "",
                         lineA_directional_traits=",".join(f"{t}:{int(v):+d}" for t, v in zip(dd["trait"], dd["direction_sign"])),
                         lineA_colocalized_traits=",".join(sorted(d["trait"])),
                         lineA_ggt_only=bool(status == "call" and traits == ["GGT"]),
                         lineA_gtex_same_sign_any=bool((dd["check_gtex_liver_same_sign_hit2"] == 1).any()),
                         lineA_not_directional_reasons="; ".join(sorted(set(d["not_directional_reason"].dropna().astype(str))))))
    genes = genes.merge(pd.DataFrame(rows), on="ensembl", how="left")
    genes["lineA_status"] = genes["lineA_status"].fillna("untestable (no row in direction table)")
    genes["lineA_sign"] = genes["lineA_sign"].fillna(0).astype(int)
    genes["lineA_call"] = genes["lineA_call"].fillna("")
    genes["published_coding_mechanism"] = genes["symbol"].map(CODING).fillna("")

    # ------------------------------------------------------------------ line B
    alt = read_backman(bdir / "GCST90083000__GCST90083000_buildGRCh38.tsv.gz")
    ast = read_backman(bdir / "GCST90083006__GCST90083006_buildGRCh38.tsv.gz")
    chen, chen_per = read_chen(zdir)
    ok = chen["chen_authors_meta_z"].notna()
    summary["chen_recomputed_vs_authors_meta"] = {
        "genes_with_ADD_row": int(len(chen)), "genes_in_authors_meta": int(ok.sum()),
        "pearson_r_z": float(np.corrcoef(chen.loc[ok, "chen_z"], chen.loc[ok, "chen_authors_meta_z"])[0, 1]),
        "share_same_sign": float((np.sign(chen.loc[ok, "chen_z"]) == np.sign(chen.loc[ok, "chen_authors_meta_z"])).mean())}

    def backman_wide(d, trait, mask):
        m = d[d["mask"] == mask].drop_duplicates("ensembl")
        return m.set_index("ensembl")[["beta", "standard_error", "p_value", "effect_allele_frequency"]].rename(
            columns=lambda c: f"{trait}_{mask}_{c}")

    for trait, d in (("ALT", alt), ("AST", ast)):
        for mask in ("M1.1", "M3.1", "M1.01"):
            genes = genes.join(backman_wide(d, trait, mask), on="ensembl")
    genes = genes.merge(chen, on="symbol", how="left")

    def designated(row, mask="M1.1"):
        t = row["driving_trait"]
        if t in ("ALT", "AST", "GGT"):
            tr = "ALT" if t == "GGT" else t
            return (f"Backman {tr} {mask}" + (" (matched trait for GGT)" if t == "GGT" else " (same trait)"),
                    row[f"{tr}_{mask}_beta"], row[f"{tr}_{mask}_p_value"])
        if t in ("NAFLD", "NASH", "PDFF"):
            return ("Chen PDFF/MASLD PTV + deleterious missense (matched trait)", row["chen_z"], row["chen_p"])
        return ("no driving trait", np.nan, np.nan)

    for mask, pre in (("M1.1", "lineB"), ("M3.1", "lineB_M3")):
        des = genes.apply(lambda r: designated(r, mask), axis=1, result_type="expand")
        genes[f"{pre}_test"], genes[f"{pre}_effect"], genes[f"{pre}_p"] = des[0], des[1], des[2]
        genes[f"{pre}_testable"] = np.isfinite(genes[f"{pre}_effect"]) & np.isfinite(genes[f"{pre}_p"])
    prim = genes["primary_evidence_class"] == "genetically_anchored"
    genes["lineB_q"] = np.nan
    genes.loc[prim, "lineB_q"] = bh(genes.loc[prim, "lineB_p"].where(genes.loc[prim, "lineB_testable"]))
    genes["lineB_M3_q_sens"] = np.nan
    genes.loc[prim, "lineB_M3_q_sens"] = bh(genes.loc[prim, "lineB_M3_p"].where(genes.loc[prim, "lineB_M3_testable"]))
    genes["lineB_q_406_sens"] = bh(genes["lineB_p"].where(genes["lineB_testable"]))
    genes["lineB_sign"] = [sign(x) if t else 0 for x, t in zip(genes["lineB_effect"], genes["lineB_testable"])]
    genes["lineB_supported"] = genes["lineB_q"] < 0.05
    genes["lineB_call"] = [call_from_burden_sign(s) if sup else "" for s, sup in zip(genes["lineB_sign"], genes["lineB_supported"])]
    genes["lineB_status"] = np.where(~genes["lineB_testable"], "untestable",
                                     np.where(genes["lineB_supported"], "supported (BH q < 0.05)", "tested, not supported"))
    genes.loc[~prim, "lineB_status"] = np.where(genes.loc[~prim, "lineB_testable"], "tested (outside primary family)", "untestable")

    # sensitivity family (i): all gene x {ALT M1.1, AST M1.1, Chen} tests of the 383 genes
    long = pd.concat([
        genes.loc[prim, ["symbol", "ensembl"]].assign(test="Backman ALT M1.1", effect=genes.loc[prim, "ALT_M1.1_beta"], p=genes.loc[prim, "ALT_M1.1_p_value"]),
        genes.loc[prim, ["symbol", "ensembl"]].assign(test="Backman AST M1.1", effect=genes.loc[prim, "AST_M1.1_beta"], p=genes.loc[prim, "AST_M1.1_p_value"]),
        genes.loc[prim, ["symbol", "ensembl"]].assign(test="Chen PDFF/MASLD", effect=genes.loc[prim, "chen_z"], p=genes.loc[prim, "chen_p"])])
    long = long[np.isfinite(long["p"])].copy()
    long["q_all_tests_sens"] = bh(long["p"])
    long.sort_values("p").to_csv(out / "lineB_all_tests_sensitivity.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ agreement
    a_ok = genes["lineA_status"] == "call"
    b_ok = genes["lineB_supported"]
    genes["agreement"] = ""
    genes.loc[a_ok & b_ok, "agreement"] = np.where(genes.loc[a_ok & b_ok, "lineA_call"] == genes.loc[a_ok & b_ok, "lineB_call"], "agree", "disagree")
    genes["headline_excluded_coding"] = (genes["published_coding_mechanism"] != "") & (genes["agreement"] == "disagree")
    genes["label"] = PROVISIONAL
    genes.sort_values(["primary_evidence_class", "symbol"]).to_csv(out / "h26_gene_direction_table.tsv", sep="\t", index=False)

    P = genes[prim]
    both_testable = (P["lineA_status"] == "call") & P["lineB_testable"]
    conc = (P.loc[both_testable, "lineA_call"] == P.loc[both_testable, "lineB_sign"].map(call_from_burden_sign))
    nom = both_testable & (P["lineB_p"] < 0.05)
    conc_nom = (P.loc[nom, "lineA_call"] == P.loc[nom, "lineB_sign"].map(call_from_burden_sign))

    def binom(k, n):
        return float(stats.binomtest(int(k), int(n), 0.5).pvalue) if n else None

    summary["primary_383"] = {
        "lineA_status": P["lineA_status"].value_counts().to_dict(),
        "lineA_calls": P.loc[P["lineA_status"] == "call", "lineA_call"].value_counts().to_dict(),
        "lineA_calls_ggt_only": int(P["lineA_ggt_only"].fillna(False).sum()),
        "lineA_calls_flagged_coding": sorted(P.loc[(P["lineA_status"] == "call") & (P["published_coding_mechanism"] != ""), "symbol"]),
        "driving_trait": P["driving_trait"].value_counts(dropna=False).to_dict(),
        "lineB_designated_test": P["lineB_test"].value_counts().to_dict(),
        "lineB_testable": int(P["lineB_testable"].sum()), "lineB_untestable": int((~P["lineB_testable"]).sum()),
        "lineB_supported_BH_q_lt_0.05": int(P["lineB_supported"].sum()),
        "lineB_supported_genes": P.loc[P["lineB_supported"], ["symbol", "lineB_test", "lineB_effect", "lineB_p", "lineB_q", "lineB_call", "lineA_status", "lineA_call"]].to_dict("records"),
        "lineB_min_q": float(P["lineB_q"].min()),
        "lineB_nominal_p_lt_0.05_count_not_significant": int((P["lineB_p"] < 0.05).sum()),
        "both_lines_supported": int(((P["lineA_status"] == "call") & P["lineB_supported"]).sum()),
        "agree": sorted(P.loc[P["agreement"] == "agree", "symbol"]),
        "disagree_headline": sorted(P.loc[(P["agreement"] == "disagree") & ~P["headline_excluded_coding"], "symbol"]),
        "disagree_coding_flagged": sorted(P.loc[(P["agreement"] == "disagree") & P["headline_excluded_coding"], "symbol"]),
        "descriptive_sign_concordance": {"n": int(both_testable.sum()), "concordant": int(conc.sum()),
                                         "binomial_p_vs_0.5": binom(conc.sum(), both_testable.sum())},
        "descriptive_sign_concordance_burden_nominal_p_lt_0.05": {
            "n": int(nom.sum()), "concordant": int(conc_nom.sum()), "binomial_p_vs_0.5": binom(conc_nom.sum(), nom.sum()),
            "note": "nominal burden p; descriptive, not a significance claim"},
        "sens_M3.1_supported": sorted(P.loc[P["lineB_M3_q_sens"] < 0.05, "symbol"]),
        "sens_all_tests_supported": sorted(set(long.loc[long["q_all_tests_sens"] < 0.05, "symbol"])),
        "sens_406_supported": sorted(genes.loc[genes["lineB_q_406_sens"] < 0.05, "symbol"]),
        "sens_excluding_ggt_only_lineA_calls": int(((P["lineA_status"] == "call") & ~P["lineA_ggt_only"].fillna(False)).sum()),
    }

    # ------------------------------------------------------------------ calibration
    cal = []
    for g, exp in CALIBRATION.items():
        ens = cat.loc[cat["symbol"] == g, "ensembl_id"]
        ens = ens.iloc[0] if len(ens) else None
        cls = cat.loc[cat["symbol"] == g, "primary_evidence_class"].iloc[0] if ens else "not in Catalog"
        row = dict(gene=g, catalog_class=cls, expected_burden_sign=exp)
        for trait, d in (("ALT", alt), ("AST", ast)):
            for mask in ("M1.1", "M3.1"):
                m = d[(d["ensembl"] == ens) & (d["mask"] == mask)]
                row[f"{trait}_{mask}_beta"] = m["beta"].iloc[0] if len(m) else np.nan
                row[f"{trait}_{mask}_p"] = m["p_value"].iloc[0] if len(m) else np.nan
        c = chen[chen["symbol"] == g]
        row["chen_z"] = c["chen_z"].iloc[0] if len(c) else np.nan
        row["chen_p"] = c["chen_p"].iloc[0] if len(c) else np.nan
        row["chen_authors_meta_z"] = c["chen_authors_meta_z"].iloc[0] if len(c) else np.nan
        la = genes[genes["symbol"] == g]
        row["lineA_status"] = la["lineA_status"].iloc[0] if len(la) else "not in the colocalized sets"
        row["lineA_call"] = la["lineA_call"].iloc[0] if len(la) else ""
        d_all = dt[dt["gene"] == g]
        row["direction_table_rows"] = ",".join(f"{t}:{int(v):+d}" for t, v in zip(d_all["trait"], d_all["direction_sign"].fillna(0)))
        b, p = row["ALT_M1.1_beta"], row["ALT_M1.1_p"]
        if exp is None:
            row["calibration"] = "recorded only"
        elif not np.isfinite(b):
            row["calibration"] = "untestable (no M1.1 ALT row)"
        elif sign(b) == exp:
            row["calibration"] = "expected sign" + (", nominal p < 0.05" if p < 0.05 else ", nominal p >= 0.05")
        else:
            row["calibration"] = "WRONG SIGN, nominal p < 0.05" if p < 0.05 else "opposite sign, nominal p >= 0.05"
        cal.append(row)
    cal = pd.DataFrame(cal)
    cal.to_csv(out / "calibration.tsv", sep="\t", index=False)
    apob = cal[cal["gene"] == "APOB"].iloc[0]
    fail = bool(cal["calibration"].str.startswith("WRONG").any() or not (apob["ALT_M1.1_beta"] > 0))
    summary["calibration"] = {"pass": not fail, "rows": cal[["gene", "catalog_class", "expected_burden_sign", "ALT_M1.1_beta", "ALT_M1.1_p", "chen_z", "chen_p", "lineA_status", "lineA_call", "calibration"]].to_dict("records")}

    # ------------------------------------------------------------------ drug targets (unverified list)
    drug = pd.DataFrame(DRUG_TARGETS, columns=["symbol", "drugs_recalled", "drug_direction", "status_recalled"])
    drug["mechanism_verified"] = False
    drug = drug.merge(cat[["symbol", "ensembl_id", "primary_evidence_class"]], on="symbol", how="left")
    drug = drug.merge(genes[["symbol", "lineA_status", "lineA_call", "lineB_status", "lineB_call", "lineB_test", "lineB_effect", "lineB_p", "lineB_q"]], on="symbol", how="left")
    drug["genetic_call_available"] = (drug["lineA_call"].fillna("") != "") | (drug["lineB_call"].fillna("") != "")
    drug["contradiction_candidate"] = [
        any(c not in ("", d) for c in (a if isinstance(a, str) else "", b if isinstance(b, str) else ""))
        for a, b, d in zip(drug["lineA_call"], drug["lineB_call"], drug["drug_direction"])]
    for trait, d in (("ALT", alt), ("AST", ast)):   # descriptive burden context for listed targets (not a call)
        drug = drug.join(backman_wide(d, trait, "M1.1"), on="ensembl_id")
    drug = drug.merge(chen[["symbol", "chen_z", "chen_p"]], on="symbol", how="left")
    drug.to_csv(out / "drug_targets_unverified.tsv", sep="\t", index=False)
    summary["drug_targets"] = {"listed": len(drug), "in_383_or_23": int(drug["lineA_status"].notna().sum()),
                               "with_supported_genetic_call": sorted(drug.loc[drug["genetic_call_available"], "symbol"]),
                               "contradiction_candidates_unverified": sorted(drug.loc[drug["contradiction_candidate"], "symbol"])}

    # ------------------------------------------------------------------ N8 (count first, direction after)
    n8 = []
    for trait, d in (("ALT", alt), ("AST", ast)):
        for mask in ("M1.1", "M3.1"):
            m = d[(d["mask"] == mask) & d["ensembl"].isin(state["ensembl_id"])].drop_duplicates("ensembl")
            n8.append(m[["ensembl", "symbol", "mask", "p_value"]].assign(trait=trait))
    n8 = pd.concat(n8)
    n8["q"] = np.nan
    for mask in ("M1.1", "M3.1"):
        k = n8["mask"] == mask
        n8.loc[k, "q"] = bh(n8.loc[k, "p_value"])
    p1 = n8[n8["mask"] == "M1.1"]
    clear = sorted(set(p1.loc[p1["q"] < 0.05, "ensembl"]))
    clear_m3 = sorted(set(n8.loc[(n8["mask"] == "M3.1") & (n8["q"] < 0.05), "ensembl"]))
    (out / "n8_count.txt").write_text(
        f"N8 ({PROVISIONAL}). Written before any burden direction was read.\n"
        f"Disease-state genes without colocalization: {len(state)}\n"
        f"Genes with a Backman M1.1 ALT or AST row: {p1['ensembl'].nunique()} ({len(p1)} gene x trait tests)\n"
        f"Genes clearing BH q < 0.05 (family = {len(p1)} tests, M1.1 pLoF, MAF <= 1%): {len(clear)}\n"
        f"Sensitivity M3.1 (own BH over {int((n8['mask'] == 'M3.1').sum())} tests): {len(clear_m3)}\n")
    betas = pd.concat([alt.assign(trait="ALT"), ast.assign(trait="AST")])[["ensembl", "mask", "trait", "beta", "standard_error", "effect_allele_frequency"]]
    n8 = n8.merge(betas, on=["ensembl", "mask", "trait"], how="left")
    n8["burden_direction_if_cleared"] = [call_from_burden_sign(sign(b)) if q < 0.05 else "" for b, q in zip(n8["beta"], n8["q"])]
    n8.sort_values(["mask", "p_value"]).to_csv(out / "n8_burden_tests.tsv", sep="\t", index=False)
    summary["N8"] = {"state_genes": len(state), "genes_tested_M1.1": int(p1["ensembl"].nunique()), "tests_M1.1": int(len(p1)),
                     "genes_clearing_BH_M1.1": len(clear),
                     "cleared_M1.1": n8[(n8["mask"] == "M1.1") & (n8["q"] < 0.05)][["symbol", "trait", "beta", "p_value", "q", "burden_direction_if_cleared"]].to_dict("records"),
                     "genes_clearing_BH_M3.1_sens": len(clear_m3),
                     "cleared_M3.1_sens": n8[(n8["mask"] == "M3.1") & (n8["q"] < 0.05)][["symbol", "trait", "beta", "p_value", "q", "burden_direction_if_cleared"]].to_dict("records")}

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating, float)):
            return None if not math.isfinite(float(o)) else float(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o

    (out / "summary.json").write_text(json.dumps(clean(summary), indent=1))
    print(json.dumps(clean(summary), indent=1))


if __name__ == "__main__":
    main()
