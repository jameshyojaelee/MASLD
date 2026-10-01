#!/usr/bin/env python3
"""A2: seal the Model A holdout (PRESPEC_A section 5; plan "Model A-final").
Reads the A1 identity crosswalk and kinship table, the per-library metadata (sex,
condition, fibrosis stage), the GSE130970 source file (controls, age at biopsy), the
GSE193066 biopsy order and the sample list of the 844-participant bulk fit.
The A1 crosswalk carries per-library identity QC built from RNA data, and this script
uses three of its columns: expr_sex (from XIST and Y-gene expression),
other_allele_fraction_at_hom (other-allele reads pooled over all homozygous identity
sites with depth >= 20; some identity sites are Model A tag positions) and
autosomal_called (number of identity sites called). It reads no per-site allele count,
no per-gene expression value and no genotype-phenotype statistic.

Cohort roles:
- GSE213621, GSE135251, GSE162694, GSE130970, GSE240729: development + 20% sealed.
- GSE174478: sealed as a whole cohort.
- GSE126848, GSE167523, PRJNA512027: development only (main effects), never ranked.
- GSE193066: excluded from every Model A allelic statistic (heavy contamination).

Library exclusions (reason codes, several may apply):
- not_in_A1_identity: in the metadata of a listed cohort but absent from A1; no individual.
- no_individual_id: A1 identity failed (too few identity sites).
- possibly_mixed: expression sex "both", or other-allele fraction at homozygous
  sites > 0.02. The fraction also contains sequencing base error, so this cut
  falls mostly on cohorts with a high base-error rate.
- sex_mismatch: metadata sex and expression sex are both M/F and differ. Metadata
  sex is empty in GSE213621, GSE135251 and GSE240729, so the check cannot act there.
- cohort_excluded_contamination: every GSE193066 library.
- stage_discordant_individual: the library belongs to an individual excluded for
  discordant stage.

Individual exclusions: home cohort GSE193066; libraries of one genetic individual
carrying different stage labels (stage_coding() from tags/t5a_allelic_table.py:
control, S0, S1, S2, stage_missing in a staged cohort, stage_not_used in a cohort
that t5a does not stage), checked over all of the individual's libraries; no library
left after the possibly_mixed and sex_mismatch exclusions.

Home cohort: the cohort of the individual's library in the earliest GEO series. The
run stops if an individual has libraries in more than one cohort and one of them is
not a GEO series (the series number of a PRJNA accession cannot be ordered against GEO).

Seal (prespecified core): within each 20% cohort, individuals left after all the
exclusions above are ranked by sha256("MASLD-A-2026-09-23" + individual_id) and the
lowest round(0.2 n) are sealed. GSE174478 individuals left after the exclusions are
all sealed. The sealed set is then closed over first-degree and same-individual
relatives (transitively), so a family is never split between development and sealed.

Frozen crosswalk (one row per metadata library of the ten cohorts): individual_id as
numbered by the A1 run given in --identity-dir, home_cohort, S / C / stage_label,
first_biopsy (True / False / unknown) and its source, autosomal_called (A1 identity
sites), run_role (development / sealed / excluded), exclusion_reason, unit_library
and unit_choice. first_biopsy: GSE193066 from the donor-keyed placement file;
GSE130970 from age_at_biopsy in the source file: True means the youngest age_at_biopsy
among this individual's GSE130970 libraries, so every individual with one GSE130970
library is True (the others False; a tie at the youngest age gives unknown); unknown
elsewhere, because no other cohort records biopsy order. unit_library is True for exactly one library per development or
sealed individual, among its libraries with run_role development or sealed, chosen
without outcomes in this order: home-cohort library, first_biopsy True > unknown >
False, most autosomal identity sites called, lowest run. unit_choice names the step
that decided.

Also writes resource_identity_defects.tsv (every same-individual or first-degree pair
with both libraries in the 844-participant bulk fit, and every GSE213621
same-individual pair with both stage labels; with the pair's IBS0 rate, and a
relation_note on first-degree pairs whose two libraries are in different cohorts) and
resource_identity_defects_by_person.tsv (every individual with more than one library
in the 844 fit). These are records, not exclusions from the 844 fit.

With --previous-seal DIR, nothing is written unless: sealed_individuals.tsv is
byte-identical to DIR's; every run -> individual_id pair in DIR's frozen crosswalk
is unchanged; every row of DIR's excluded_libraries.tsv is present with the same
reason and every added row is not_in_A1_identity; and every per-cohort individual,
role and stage count in DIR's a2_summary.json is unchanged. Whether
frozen_crosswalk.tsv and excluded_libraries.tsv are byte-identical to DIR's is
recorded, not required.

Rules for later steps (not enforced here): join through frozen_crosswalk.tsv, checked
against frozen_crosswalk.sha256; never take individual IDs from a newer identity run.
The development set is the crosswalk rows with run_role development and unit_library
True, and no run listed in excluded_libraries.tsv. Stage eligibility follows
unit_library; first_biopsy unknown is eligible.

The run writes a byte copy of this script and its sha256 into the output directory.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd

SALT = "MASLD-A-2026-09-23"
SEAL_FRACTION = 0.2
SEAL_20PCT = ["GSE213621", "GSE135251", "GSE162694", "GSE130970", "GSE240729"]
SEAL_WHOLE = ["GSE174478"]
DEV_ONLY = ["GSE126848", "GSE167523", "PRJNA512027"]
EXCLUDED_COHORTS = {"GSE193066": "cohort_excluded_contamination"}
ALL_COHORTS = SEAL_20PCT + SEAL_WHOLE + DEV_ONLY + list(EXCLUDED_COHORTS)
MIXED_OTHER_ALLELE_FRACTION = 0.02
FAMILY = ["first_degree", "same_individual"]
LIB_CODES = ["not_in_A1_identity", "no_individual_id", "possibly_mixed", "sex_mismatch",
             "cohort_excluded_contamination", "stage_discordant_individual"]
IND_CODES = ["cohort_excluded_contamination", "stage_discordant_same_individual",
             "no_usable_library:possibly_mixed", "no_usable_library:possibly_mixed+sex_mismatch",
             "no_usable_library:sex_mismatch"]
STAGE_KEYS = ["S0", "S1", "S2", "control", "stage_missing", "stage_not_used"]
T5A = Path(__file__).resolve().parent / "tags" / "t5a_allelic_table.py"
FIRST_RANK = {"True": 0, "unknown": 1, "False": 2}
OUT_FILES = ["sealed_individuals.tsv", "excluded_libraries.tsv", "frozen_crosswalk.tsv",
             "frozen_crosswalk.sha256", "resource_identity_defects.tsv",
             "resource_identity_defects_by_person.tsv", "a2_summary.json", "a2_seal_holdout.py"]


def h(x):
    return hashlib.sha256((SALT + x).encode()).hexdigest()


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(p):
    return sha256_bytes(Path(p).read_bytes())


def load_t5a():
    spec = importlib.util.spec_from_file_location("t5a_allelic_table", T5A)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def stage_label(s, c, staged):
    if c == 1:
        return "control"
    if pd.notna(s):
        return f"S{int(s)}"
    return "stage_missing" if staged else "stage_not_used"


def join_codes(codes, order):
    if set(codes) - set(order):
        raise SystemExit(f"unlisted reason codes: {sorted(set(codes) - set(order))}")
    return ";".join(k for k in order if k in codes)


def count_codes(values):
    out = {}
    for v in values:
        for k in (v.split(";") if v else []):
            out[k] = out.get(k, 0) + 1
    return out


def check_previous(prev, files, xw, per):
    """Compare this run with an earlier seal directory. Returns the check record;
    stops the run on any difference."""
    prev = Path(prev)
    fails = []
    old_sha = sha256_file(prev / "sealed_individuals.tsv")
    new_sha = sha256_bytes(files["sealed_individuals.tsv"].encode())
    if old_sha != new_sha:
        fails.append(f"sealed_individuals.tsv differs: {old_sha} vs {new_sha}")
    stated = (prev / "frozen_crosswalk.sha256").read_text().split()[0]
    old_xw_sha = sha256_file(prev / "frozen_crosswalk.tsv")
    if stated != old_xw_sha:
        fails.append("earlier frozen_crosswalk.tsv does not match its frozen_crosswalk.sha256")
    old_xw = pd.read_csv(prev / "frozen_crosswalk.tsv", sep="\t", dtype=str, keep_default_na=False)
    new_id = xw.set_index("run")["individual_id"].fillna("").astype(str)
    missing = sorted(set(old_xw["run"]) - set(new_id.index))
    changed = sorted(r for r, i in zip(old_xw["run"], old_xw["individual_id"]) if r in new_id.index and new_id[r] != i)
    added = sorted(set(new_id.index) - set(old_xw["run"]))
    reason = xw.set_index("run")["exclusion_reason"]
    added_bad = [r for r in added if reason[r] != "not_in_A1_identity"]
    if missing or changed or added_bad:
        fails.append(f"crosswalk runs missing {missing[:5]}, individual_id changed {changed[:5]}, "
                     f"added without not_in_A1_identity {added_bad[:5]}")
    old_ex = pd.read_csv(prev / "excluded_libraries.tsv", sep="\t", dtype=str, keep_default_na=False)
    ex_changed = sorted(r for r, k in zip(old_ex["run"], old_ex["reason"]) if reason.get(r) != k)
    if ex_changed:
        fails.append(f"excluded_libraries reasons changed: {ex_changed[:5]}")
    old_per = json.loads((prev / "a2_summary.json").read_text())["per_cohort"]
    cmp_fail = []
    for c, o in old_per.items():
        n = per[c]
        for k in ["individuals", "development", "sealed", "excluded"]:
            if o[k] != n[k]:
                cmp_fail.append(f"{c} {k}")
        if o["excluded_by_reason"] != n["excluded_by_reason"]:
            cmp_fail.append(f"{c} excluded_by_reason")
        for role in ["development_by_stage", "sealed_by_stage"]:
            ns = dict(n[role])
            if "unstaged" in o[role]:  # count format of the 2026-09-26 seal
                ns["unstaged"] = ns.pop("stage_missing") + ns.pop("stage_not_used")
            if o[role] != ns:
                cmp_fail.append(f"{c} {role}")
    if cmp_fail:
        fails.append(f"per-cohort counts differ: {cmp_fail}")
    if fails:
        raise SystemExit("previous-seal check failed; nothing written:\n" + "\n".join(fails))
    return {"previous_seal_dir": str(prev.resolve()),
            "sealed_individuals_sha256_previous": old_sha, "sealed_individuals_sha256_now": new_sha,
            "sealed_individuals_byte_identical": True,
            "previous_frozen_crosswalk_sha256": old_xw_sha,
            "frozen_crosswalk_byte_identical": old_xw_sha == sha256_bytes(files["frozen_crosswalk.tsv"].encode()),
            "excluded_libraries_byte_identical":
                sha256_file(prev / "excluded_libraries.tsv") == sha256_bytes(files["excluded_libraries.tsv"].encode()),
            "runs_in_previous_crosswalk": int(len(old_xw)), "run_to_individual_id_unchanged": True,
            "runs_added": added, "previous_excluded_libraries_unchanged": int(len(old_ex)),
            "per_cohort_role_and_stage_counts_unchanged": True}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--identity-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--gse130970-controls", required=True, help="GSE130970_source_diagnosis.tsv")
    ap.add_argument("--gse193066-placement", required=True)
    ap.add_argument("--fit-samples", required=True, help="model_design.tsv of the 844-participant bulk fit")
    ap.add_argument("--previous-seal", help="earlier seal directory that this run must reproduce")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    script_bytes = Path(__file__).resolve().read_bytes()
    out = Path(a.out)
    present = [f for f in OUT_FILES if (out / f).exists()]
    if present:
        raise SystemExit(f"refusing to overwrite existing seal outputs in {out}: {present}")

    idt = Path(a.identity_dir)
    x = pd.read_csv(idt / "sample_to_individual.tsv", sep="\t")
    kin = pd.read_csv(idt / "kinship_related_pairs.tsv", sep="\t")
    meta = pd.read_csv(a.metadata, keep_default_na=False, na_values=[""])
    c130 = pd.read_csv(a.gse130970_controls, sep="\t")
    placement = pd.read_csv(a.gse193066_placement, sep="\t")
    fit = set(pd.read_csv(a.fit_samples, sep="\t")["sample_id"])
    if not x["run"].is_unique or meta["sample_id"].duplicated().any():
        raise SystemExit("duplicate runs in sample_to_individual.tsv or the metadata")
    n_lib = len(x)
    unknown = (set(x["cohort"]) | set(meta["dataset"])) - set(ALL_COHORTS)
    if unknown:
        raise SystemExit(f"cohorts without a role: {sorted(unknown)}")

    # harmonized stage per library (t5a stage_coding: GSE213621 from its condition label,
    # GSE130970 controls = the 6 source controls)
    t5a = load_t5a()
    stage = t5a.stage_coding(meta, set(c130.loc[c130["source_control_status"] == "Control", "run"]))
    stage = stage.merge(meta[["sample_id", "dataset", "condition", "fibrosis_stage"]]
                        .rename(columns={"sample_id": "run", "dataset": "cohort"}), on="run", how="left")
    stage["stage_label"] = [stage_label(s, c, coh in t5a.STAGED)
                            for s, c, coh in zip(stage["S"], stage["C"], stage["cohort"])]
    x = x.merge(stage[["run", "S", "C", "stage_label", "condition", "fibrosis_stage"]], on="run", how="left")
    if x["C"].isna().any():
        raise SystemExit(f"libraries without metadata: {x.loc[x['C'].isna(), 'run'].tolist()}")
    if len(x) != n_lib:
        raise SystemExit("metadata join changed the number of libraries")
    if (x["run"].map(stage.set_index("run")["cohort"]) != x["cohort"]).any():
        raise SystemExit("A1 cohort and metadata dataset disagree for some runs")
    not_in_a1 =stage[~stage["run"].isin(x["run"])].copy()

    # library-level exclusions
    codes = {r: set() for r in x["run"]}
    has_id = x["individual_id"].notna()
    mixed = (x["expr_sex"] == "both") | (x["other_allele_fraction_at_hom"] > MIXED_OTHER_ALLELE_FRACTION)
    sexmm = x["metadata_sex"].isin(["M", "F"]) & x["expr_sex"].isin(["M", "F"]) & (x["metadata_sex"] != x["expr_sex"])
    for r in x.loc[~has_id, "run"]:
        codes[r].add("no_individual_id")
    for r in x.loc[mixed, "run"]:
        codes[r].add("possibly_mixed")
    for r in x.loc[sexmm, "run"]:
        codes[r].add("sex_mismatch")
    for coh, code in EXCLUDED_COHORTS.items():
        for r in x.loc[x["cohort"] == coh, "run"]:
            codes[r].add(code)

    # individuals: home cohort = cohort of the earliest-series library
    xi = x[has_id]
    cohorts_of = xi.groupby("individual_id")["cohort"].agg(lambda s: sorted(set(s)))
    multi = cohorts_of[cohorts_of.map(len) > 1]
    non_geo = multi[multi.map(lambda cs: any(not c.startswith("GSE") for c in cs))]
    if len(non_geo):
        raise SystemExit(f"individuals in more than one cohort, one not a GEO series: {non_geo.to_dict()}")
    n_multi_cohort = int(len(multi))
    home = xi.assign(num=xi["cohort"].str.extract(r"(\d+)")[0].astype(int)).sort_values(["num", "run"]) \
             .drop_duplicates("individual_id").set_index("individual_id")["cohort"]
    ind_codes = {i: set() for i in home.index}
    for i, c in home.items():
        if c in EXCLUDED_COHORTS:
            ind_codes[i].add(EXCLUDED_COHORTS[c])
    labels = xi.groupby("individual_id")["stage_label"].agg(lambda s: sorted(set(s)))
    for i in labels.index[labels.map(len) > 1]:
        ind_codes[i].add("stage_discordant_same_individual")
        for r in xi.loc[xi["individual_id"] == i, "run"]:
            codes[r].add("stage_discordant_individual")
    lib_level = {"possibly_mixed", "sex_mismatch"}
    for i, g in xi.groupby("individual_id"):
        runs = g["run"].tolist()
        if all(codes[r] & lib_level for r in runs):
            ind_codes[i].add("no_usable_library:" + "+".join(sorted(set().union(*(codes[r] & lib_level for r in runs)))))
    excluded = {i for i, k in ind_codes.items() if k}
    eligible = home[~home.index.isin(excluded)]

    # seal: lowest 20% by salted hash per cohort; GSE174478 whole
    role, reason, rank_info = {}, {}, {}
    for c in SEAL_20PCT:
        inds = sorted(eligible[eligible == c].index, key=h)
        n_seal = round(SEAL_FRACTION * len(inds))
        rank_info[c] = {"ranked_n": len(inds), "seal_fraction_times_n": round(SEAL_FRACTION * len(inds), 4),
                        "n_seal_hash": n_seal}
        for i in inds[:n_seal]:
            role[i], reason[i] = "sealed", "hash_lowest_20pct"
    for c in SEAL_WHOLE:
        for i in eligible[eligible == c].index:
            role[i], reason[i] = "sealed", "whole_cohort_sealed"

    # close the sealed set over first-degree and same-individual relatives (transitive,
    # through any individual, so a family is never split)
    run_ind = x.set_index("run")["individual_id"]
    fam = {}
    rel = kin[kin["relation"].isin(FAMILY)]
    for a_, b_ in zip(rel["#IID1"].map(run_ind), rel["IID2"].map(run_ind)):
        if pd.notna(a_) and pd.notna(b_) and a_ != b_:
            fam.setdefault(a_, set()).add(b_); fam.setdefault(b_, set()).add(a_)
    reached = {i for i, r in role.items() if r == "sealed"}
    frontier = list(reached)
    while frontier:
        i = frontier.pop()
        for j in fam.get(i, ()):
            if j not in reached:
                reached.add(j); frontier.append(j)
    moved = sorted(j for j in reached if j in eligible.index and role.get(j) != "sealed")
    for j in moved:
        role[j], reason[j] = "sealed", "relative_of_sealed"

    for i, c in eligible.items():
        if i not in role:
            role[i] = "development"
            reason[i] = "hash_rank_above_20pct" if c in SEAL_20PCT else "development_only_cohort"
    for i in excluded:
        role[i], reason[i] = "excluded", join_codes(ind_codes[i], IND_CODES)

    tab = pd.DataFrame({"individual_id": sorted(home.index)})
    tab["home_cohort"] = tab["individual_id"].map(home)
    tab["role"] = tab["individual_id"].map(role)
    tab["reason"] = tab["individual_id"].map(reason)
    if tab["role"].isna().any():
        raise SystemExit("individual without a role")

    # moved relatives: each kinship row that links a moved individual to a sealed one
    moved_rec = []
    for r1, r2, rl, k_, ibs0, nsnp in zip(rel["#IID1"], rel["IID2"], rel["relation"], rel["KINSHIP"],
                                          rel["IBS0"], rel["NSNP"]):
        i1, i2 = run_ind.get(r1), run_ind.get(r2)
        for (mi, mr), (si, sr) in (((i1, r1), (i2, r2)), ((i2, r2), (i1, r1))):
            if mi in moved and si is not None and pd.notna(si) and role.get(si) == "sealed" and si != mi:
                moved_rec.append({"individual_id": mi, "home_cohort": home[mi], "run": mr,
                                  "sealed_relative": si, "relative_home_cohort": home[si],
                                  "relative_reason": reason[si], "relative_run": sr, "relation": rl,
                                  "kinship": float(k_), "ibs0": int(ibs0), "nsnp": int(nsnp),
                                  "ibs0_rate": round(int(ibs0) / int(nsnp), 4)})

    # per-library role, first biopsy, unit library
    x["exclusion_reason"] = x["run"].map(lambda r: join_codes(codes[r], LIB_CODES))
    ind_role = tab.set_index("individual_id")["role"]
    x["home_cohort"] = x["individual_id"].map(home)
    x["run_role"] = x["individual_id"].map(ind_role).where(x["exclusion_reason"] == "", "excluded")
    if x["run_role"].isna().any():
        raise SystemExit("library without a role")
    orphan = x[(x["individual_id"].map(ind_role) == "excluded") & (x["exclusion_reason"] == "")]
    if len(orphan):
        raise SystemExit(f"libraries of excluded individuals without a library reason: {orphan['run'].tolist()}")

    x["first_biopsy"], x["first_biopsy_source"] = "unknown", "not_recorded"
    first = placement.set_index("sample_id")["is_first_biopsy"].astype(bool)
    g193 = (x["cohort"] == "GSE193066") & x["run"].isin(first.index)
    x.loc[g193, "first_biopsy"] = x.loc[g193, "run"].map(first).map({True: "True", False: "False"})
    x.loc[g193, "first_biopsy_source"] = "gse193066_donor_keyed_placement"
    age = c130.set_index("run")["age_at_biopsy"]
    g130 = x[(x["cohort"] == "GSE130970") & x["individual_id"].notna()]
    for i, g in g130.groupby("individual_id"):
        ages = g["run"].map(age)
        if ages.isna().any():
            continue
        youngest = ages == ages.min()
        val = ["unknown" if (y and youngest.sum() > 1) else ("True" if y else "False") for y in youngest]
        x.loc[g.index, "first_biopsy"] = val
        x.loc[g.index, "first_biopsy_source"] = "GSE130970_source_age_at_biopsy"

    usable = x[x["run_role"].isin(["development", "sealed"])].copy()
    usable["is_home"] = usable["cohort"] == usable["home_cohort"]
    usable["fb_rank"] = usable["first_biopsy"].map(FIRST_RANK)
    usable = usable.sort_values(["individual_id", "is_home", "fb_rank", "autosomal_called", "run"],
                                ascending=[True, False, True, False, True])
    choice = {}
    for i, g in usable.groupby("individual_id", sort=False):
        top = g.iloc[0]
        if len(g) == 1:
            step = "only_usable_library"
        else:
            nxt = g.iloc[1]
            step = ("home_cohort" if top["is_home"] != nxt["is_home"] else
                    "first_biopsy" if top["fb_rank"] != nxt["fb_rank"] else
                    "most_identity_sites_called" if top["autosomal_called"] != nxt["autosomal_called"] else
                    "lowest_run")
        choice[top["run"]] = step
    x["unit_library"] = x["run"].isin(choice)
    x["unit_choice"] = x["run"].map(choice).fillna("")
    kept = set(tab.loc[tab["role"].isin(["development", "sealed"]), "individual_id"])
    n_unit = x[x["unit_library"]].groupby("individual_id").size()
    if set(n_unit.index) != kept or (n_unit != 1).any():
        raise SystemExit("unit_library is not exactly one per development or sealed individual")

    # the metadata libraries that A1 never saw
    extra = not_in_a1.assign(individual_id=pd.NA, home_cohort=pd.NA, autosomal_called=pd.NA,
                             run_role="excluded", exclusion_reason="not_in_A1_identity",
                             first_biopsy="unknown", first_biopsy_source="not_recorded",
                             unit_library=False, unit_choice="")
    xw_cols = ["run", "cohort", "individual_id", "home_cohort", "S", "C", "stage_label", "first_biopsy",
               "first_biopsy_source", "autosomal_called", "run_role", "exclusion_reason", "unit_library",
               "unit_choice"]
    xw = pd.concat([x[xw_cols], extra[xw_cols]], ignore_index=True).sort_values(["cohort", "run"])
    xw["S"] = xw["S"].astype("Int64")
    xw["C"] = xw["C"].astype(int)
    xw["autosomal_called"] = xw["autosomal_called"].astype("Int64")
    exl = xw.loc[xw["exclusion_reason"] != "", ["run", "cohort", "exclusion_reason"]] \
            .rename(columns={"exclusion_reason": "reason"})

    # identity records for the Resource (not exclusions from the 844 fit)
    lab, cond, coh = (x.set_index("run")[k] for k in ("stage_label", "condition", "cohort"))
    d = rel.rename(columns={"#IID1": "run_1", "IID2": "run_2", "KINSHIP": "kinship", "IBS0": "ibs0",
                            "NSNP": "nsnp"})[["run_1", "run_2", "relation", "kinship", "ibs0", "nsnp"]].copy()
    d["ibs0_rate"] = (d["ibs0"] / d["nsnp"]).round(4)
    for k in ("1", "2"):
        d[f"cohort_{k}"] = d[f"run_{k}"].map(coh)
        d[f"individual_id_{k}"] = d[f"run_{k}"].map(run_ind)
        d[f"condition_{k}"] = d[f"run_{k}"].map(cond)
        d[f"stage_{k}"] = d[f"run_{k}"].map(lab)
        d[f"in_844_fit_{k}"] = d[f"run_{k}"].isin(fit)
    d["stage_concordant"] = d["stage_1"] == d["stage_2"]
    # a first-degree call between libraries of two cohorts is not a documented family:
    # note it with the IBS0 rates of the other related pairs outside the excluded cohorts
    kr = kin.assign(c1=kin["#IID1"].map(coh), c2=kin["IID2"].map(coh), rate=kin["IBS0"] / kin["NSNP"])
    clean = ~kr["c1"].isin(list(EXCLUDED_COHORTS)) & ~kr["c2"].isin(list(EXCLUDED_COHORTS))
    ref1 = kr.loc[clean & (kr["relation"] == "first_degree") & (kr["c1"] == kr["c2"]), "rate"]
    ref2 = kr.loc[clean & (kr["relation"] == "second_degree"), "rate"]
    outside = "outside " + ", ".join(EXCLUDED_COHORTS)
    ibs0_reference = {"first_degree_within_one_cohort_" + outside.replace(" ", "_"):
                          {"n_pairs": int(len(ref1)), "min": round(float(ref1.min()), 4), "max": round(float(ref1.max()), 4)},
                      "second_degree_" + outside.replace(" ", "_"):
                          {"n_pairs": int(len(ref2)), "min": round(float(ref2.min()), 4), "max": round(float(ref2.max()), 4)}}
    cross1 = (d["relation"] == "first_degree") & (d["cohort_1"] != d["cohort_2"])
    d["relation_note"] = ""
    d.loc[cross1, "relation_note"] = [
        f"kinship in first-degree range, across two cohorts, not confirmed; IBS0 rate {r:.4f} "
        f"(first-degree pairs within one cohort {outside}: {ref1.min():.4f}-{ref1.max():.4f}, n={len(ref1)}; "
        f"second-degree pairs {outside}: {ref2.min():.4f}-{ref2.max():.4f}, n={len(ref2)})"
        for r in d.loc[cross1, "ibs0_rate"]]
    both_fit = d["in_844_fit_1"] & d["in_844_fit_2"]
    g213 = (d["relation"] == "same_individual") & (d["cohort_1"] == "GSE213621") & (d["cohort_2"] == "GSE213621")
    d["record"] = [";".join(k for k, f in (("both_in_844_fit", bf), ("GSE213621_same_individual", gf)) if f)
                   for bf, gf in zip(both_fit, g213)]
    d = d[both_fit | g213].sort_values(["relation", "cohort_1", "run_1"])
    d = d[["record", "relation", "relation_note", "kinship", "ibs0", "nsnp", "ibs0_rate", "run_1", "run_2",
           "cohort_1", "cohort_2", "individual_id_1",
           "individual_id_2", "in_844_fit_1", "in_844_fit_2", "condition_1", "condition_2", "stage_1",
           "stage_2", "stage_concordant"]]

    # people in the 844-participant fit
    xf = x[x["run"].isin(fit)]
    xfi = xf[xf["individual_id"].notna()]
    per_ind = xfi.groupby("individual_id")
    dup = per_ind.filter(lambda g: len(g) > 1)
    dp = dup.groupby("individual_id").agg(
        cohort=("cohort", lambda s: ";".join(sorted(set(s)))),
        n_libraries_in_844_fit=("run", "size"),
        runs=("run", lambda s: ";".join(sorted(s))),
        conditions=("condition", lambda s: ";".join(str(v) for v in s)),
        stage_labels=("stage_label", lambda s: ";".join(s)),
        stage_concordant=("stage_label", lambda s: s.nunique() == 1)).reset_index()
    sizes = per_ind.size()
    fit_people = {
        "fit_libraries": len(fit),
        "fit_libraries_in_A1": int(len(xf)),
        "fit_libraries_not_in_A1": int(len(fit - set(x["run"]))),
        "fit_libraries_without_individual_id": int(xf["individual_id"].isna().sum()),
        "fit_libraries_without_individual_id_by_cohort": xf[xf["individual_id"].isna()]["cohort"].value_counts().to_dict(),
        "identified_individuals": int(len(sizes)),
        "individuals_with_more_than_one_library": int((sizes > 1).sum()),
        "individuals_by_library_count": {str(k): int(v) for k, v in sizes[sizes > 1].value_counts().sort_index().items()},
        "surplus_libraries": int((sizes - 1).sum()),
        "people_at_most": int(len(sizes) + xf["individual_id"].isna().sum() + len(fit - set(x["run"]))),
        "stage_discordant_libraries_in_fit": int(xf["exclusion_reason"].str.contains("stage_discordant_individual").sum()),
        "individuals_with_discordant_stage_in_fit": int((~dp["stage_concordant"]).sum()),
    }

    # summary (every table is built before the first write)
    ind_lab = xi.groupby("individual_id")["stage_label"].first()
    unit_row = x[x["unit_library"]].set_index("individual_id")
    tl = tab.assign(stage_label=tab["individual_id"].map(ind_lab))
    per = {}
    for c in ALL_COHORTS:
        t = tl[tl["home_cohort"] == c]
        ex = t[t["role"] == "excluded"]
        lc = x[x["cohort"] == c]
        uc = unit_row[unit_row["home_cohort"] == c]
        per[c] = {
            "libraries_in_metadata": int((meta["dataset"] == c).sum()),
            "libraries_in_A1": int(len(lc)),
            "libraries_excluded": int((xw[xw["cohort"] == c]["exclusion_reason"] != "").sum()),
            "libraries_excluded_by_reason": count_codes(xw[xw["cohort"] == c]["exclusion_reason"]),
            "libraries_with_metadata_and_expression_sex_MF": int(
                (lc["metadata_sex"].isin(["M", "F"]) & lc["expr_sex"].isin(["M", "F"])).sum()),
            "individuals": int(len(t)),
            "development": int((t["role"] == "development").sum()),
            "sealed": int((t["role"] == "sealed").sum()),
            "excluded": int(len(ex)),
            "excluded_by_reason": count_codes(ex["reason"]),
            "development_by_stage": {k: int(((t["role"] == "development") & (t["stage_label"] == k)).sum())
                                     for k in STAGE_KEYS},
            "sealed_by_stage": {k: int(((t["role"] == "sealed") & (t["stage_label"] == k)).sum()) for k in STAGE_KEYS},
            "sealed_by_reason": t[t["role"] == "sealed"]["reason"].value_counts().to_dict(),
            "moved_relatives": int(sum(1 for j in moved if home[j] == c)),
            "unit_libraries": int(len(uc)),
            "unit_choice": uc["unit_choice"].value_counts().to_dict(),
        }
        if c in rank_info:
            per[c].update(rank_info[c])
        if c in DEV_ONLY:
            dev_units = uc[uc["run_role"] == "development"]
            per[c]["development_recorded_condition"] = dev_units["condition"].fillna("NA").astype(str) \
                .value_counts().to_dict()
            per[c]["development_recorded_fibrosis_stage"] = dev_units["fibrosis_stage"] \
                .map(lambda v: "NA" if pd.isna(v) else str(int(v))).value_counts().to_dict()

    # libraries coded as control (C = 1) whose metadata still records a fibrosis stage above 0
    ctl = xw.merge(stage[["run", "condition", "fibrosis_stage"]], on="run", how="left")
    ctl = ctl[(ctl["C"] == 1) & (ctl["fibrosis_stage"] > 0)]
    nas = meta.set_index("sample_id")["nas_score"]
    controls_with_stage = [
        {"run": r.run, "cohort": r.cohort,
         "individual_id": None if pd.isna(r.individual_id) else str(r.individual_id),
         "condition": str(r.condition), "fibrosis_stage": float(r.fibrosis_stage),
         "nas_score": None if pd.isna(nas.get(r.run)) else float(nas.get(r.run)),
         "stage_label": r.stage_label, "run_role": r.run_role, "unit_library": bool(r.unit_library)}
        for r in ctl.itertuples()]

    files = {
        "sealed_individuals.tsv": tab.to_csv(sep="\t", index=False),
        "excluded_libraries.tsv": exl.to_csv(sep="\t", index=False),
        "frozen_crosswalk.tsv": xw.to_csv(sep="\t", index=False),
        "resource_identity_defects.tsv": d.to_csv(sep="\t", index=False),
        "resource_identity_defects_by_person.tsv": dp.to_csv(sep="\t", index=False),
    }
    prev_check = check_previous(a.previous_seal, files, xw, per) if a.previous_seal else None

    summ = {
        "source_script": str(Path(__file__).resolve()),
        "source_script_sha256": sha256_bytes(script_bytes),
        "source_script_bytes": len(script_bytes),
        "inputs": {str(p): sha256_file(p) for p in [idt / "sample_to_individual.tsv", idt / "kinship_related_pairs.tsv",
                                                    Path(a.metadata), Path(a.gse130970_controls),
                                                    Path(a.gse193066_placement), Path(a.fit_samples), T5A]},
        "rules": {"salt": SALT, "seal_fraction": SEAL_FRACTION, "n_seal": "round(seal_fraction * ranked_n)",
                  "seal_20pct_cohorts": SEAL_20PCT, "seal_whole_cohorts": SEAL_WHOLE,
                  "development_only_cohorts": DEV_ONLY, "excluded_cohorts": EXCLUDED_COHORTS,
                  "possibly_mixed": f"expr_sex == 'both' or other_allele_fraction_at_hom > {MIXED_OTHER_ALLELE_FRACTION}",
                  "sex_mismatch": "metadata_sex and expr_sex both in {M, F} and differ",
                  "stage_label": "stage_coding() in " + str(T5A) + "; control if C == 1, else S0/S1/S2; "
                                 "stage_missing if the cohort is in t5a STAGED, else stage_not_used",
                  "stage_discordance": "over all of an individual's libraries",
                  "ranking_population": "individuals not excluded, by home cohort (after every exclusion)",
                  "relative_closure": f"relatives by {FAMILY} of a sealed individual are sealed (transitive)",
                  "first_biopsy": "GSE193066 placement file; GSE130970 True = youngest age_at_biopsy among the "
                                  "individual's GSE130970 libraries (every one-library individual is True); "
                                  "else unknown (biopsy order not recorded)",
                  "stage_eligibility_for_later_steps": "follows unit_library; first_biopsy unknown is eligible",
                  "unit_library": "per development or sealed individual among usable libraries: home cohort, "
                                  "first_biopsy True > unknown > False, most autosomal_called, lowest run"},
        "totals": {"libraries_in_metadata": int(len(meta)), "libraries_in_A1": int(len(x)),
                   "crosswalk_rows": int(len(xw)), "libraries_excluded": int(len(exl)),
                   "libraries_excluded_by_reason": count_codes(exl["reason"]),
                   "individuals": int(len(tab)),
                   "development": int((tab["role"] == "development").sum()),
                   "sealed": int((tab["role"] == "sealed").sum()),
                   "excluded": int((tab["role"] == "excluded").sum()),
                   "individuals_in_more_than_one_cohort": n_multi_cohort,
                   "moved_into_sealed_as_relatives": len(moved),
                   "sealed_by_reason": tab[tab["role"] == "sealed"]["reason"].value_counts().to_dict(),
                   "unit_libraries": int(x["unit_library"].sum()),
                   "individuals_with_more_than_one_usable_library":
                       int((usable.groupby("individual_id").size() > 1).sum())},
        "moved_relatives": moved_rec,
        "per_cohort": per,
        "resource_identity_defects": {
            "rows": int(len(d)),
            "both_in_844_fit": int(d["record"].str.contains("both_in_844_fit").sum()),
            "both_in_844_fit_by_relation": d[d["record"].str.contains("both_in_844_fit")]["relation"]
            .value_counts().to_dict(),
            "GSE213621_same_individual": int(d["record"].str.contains("GSE213621_same_individual").sum()),
            "GSE213621_same_individual_stage_discordant": int(
                (d["record"].str.contains("GSE213621_same_individual") & ~d["stage_concordant"]).sum()),
            "ibs0_rate_reference": ibs0_reference,
            "first_degree_across_two_cohorts": d.loc[d["relation_note"] != "",
                                                     ["run_1", "run_2", "cohort_1", "cohort_2", "kinship",
                                                      "ibs0", "nsnp", "ibs0_rate", "relation_note"]]
            .to_dict(orient="records"),
            "people_in_844_fit": fit_people},
        "controls_with_recorded_fibrosis_stage_above_0": {
            "rule": "coded as control (C = 1, stage_label control), not staged; used only for main effects",
            "libraries": controls_with_stage},
        "previous_seal_check": prev_check,
    }

    out.mkdir(parents=True, exist_ok=True)
    for f, s in files.items():
        (out / f).write_text(s)
    (out / "a2_seal_holdout.py").write_bytes(script_bytes)
    sha_xw = sha256_file(out / "frozen_crosswalk.tsv")
    (out / "frozen_crosswalk.sha256").write_text(f"{sha_xw}  frozen_crosswalk.tsv\n")
    summ["outputs"] = {f: sha256_file(out / f) for f in list(files) + ["a2_seal_holdout.py"]}
    if summ["outputs"]["a2_seal_holdout.py"] != summ["source_script_sha256"]:
        raise SystemExit("script copy does not match the script that ran")
    summ["sha256_sealed_individuals"] = summ["outputs"]["sealed_individuals.tsv"]
    summ["sha256_frozen_crosswalk"] = sha_xw
    (out / "a2_summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
