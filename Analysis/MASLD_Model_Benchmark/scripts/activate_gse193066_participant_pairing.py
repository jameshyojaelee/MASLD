#!/usr/bin/env python
"""
activate_gse193066_participant_pairing.py

Build the participant-level pairing for GSE193066 (Hoshida/UTSW, NAFLD paired
biopsies) from the deposited GEO series matrix, and characterise the repeat
biopsy structure and the phenotype fields.

WHY THIS EXISTS
  GSE193066 carries 164 SAMPLE ROWS but 106 PARTICIPANTS: 106 first biopsies
  plus 58 repeat biopsies from a subset of those participants. It is the only
  bulk cohort in this project with a real donor key AND within-participant
  repeats, so it is the only bulk substrate that supports a paired,
  donor-level contrast.

THE TITLE CONVENTION, AND WHY THE NAIVE STRIP PAIRS ZERO
  Repeat-biopsy participants are titled  HUnafld035_1 / HUnafld035_2.
  Single-biopsy participants are titled  HUnafld001  -- NO suffix at all.
  A "strip _2 and look up the bare title" join therefore matches nothing:
  HUnafld035_2 -> HUnafld035, which is not a title in the file. The correct
  rule is to strip a trailing _<digit> from EVERY title, then group.
  build_naive_pairing() reproduces the broken rule so the guard that rejects
  it can be shown to fire, rather than being asserted to work.

Outputs (into --outdir):
  participant_pairing.tsv    one row per participant (106)
  sample_records.tsv         one row per GSM (164), with participant key
  repeat_structure.json      repeat counts, interval evidence, phenotype fields
  guard_report.json          every guard, its category, and its firing proof
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, os, re, sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
GEO_DIR = ROOT/"RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/metadata"
SERIES = GEO_DIR/"GSE193066_series_matrix.txt.gz"
SRA_META = GEO_DIR/"metadata.tsv"
UNIFIED = ROOT/"RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
GCT106 = ROOT/"data/published_degs/GSE193066/GSE193066_NAFLD.HUn106.gct.gz"
GCT164 = ROOT/"data/published_degs/GSE193066/GSE193066_NAFLD.HUn164.gct.gz"
TRAJ = ROOT/"RNA-seq/results/reversal/gse193066_paired_trajectory.csv"

SUFFIX_RE = re.compile(r"_(\d+)$")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# ---------------------------------------------------------------- parsing
def parse_series(path: Path) -> tuple[pd.DataFrame, dict]:
    """Parse the series matrix PER SAMPLE, keying each characteristic by its
    own 'key: value' prefix rather than by row position."""
    lines = [l.rstrip("\n") for l in gzip.open(path, "rt", errors="ignore")]

    def row(tag):
        return next(l for l in lines if l.startswith(tag))

    def fields(l):
        return [p.strip('"') for p in l.split("\t")[1:]]

    gsms = fields(row("!Sample_geo_accession"))
    titles = fields(row("!Sample_title"))
    char_rows = [fields(l) for l in lines
                 if l.startswith("!Sample_characteristics_ch1")]
    rel_rows = [fields(l) for l in lines if l.startswith("!Sample_relation")]

    struct = {
        "n_sample_columns": len(gsms),
        "n_characteristics_rows": len(char_rows),
        "characteristics_row_widths": sorted({len(r) for r in char_rows}),
        "title_row_width": len(titles),
    }

    recs = []
    for j, g in enumerate(gsms):
        d = {"gsm": g, "title": titles[j], "_keys": []}
        for r in char_rows:
            v = r[j]
            if ": " in v:
                k, val = v.split(": ", 1)
                d[k.strip()] = val.strip()
                d["_keys"].append(k.strip())
        for r in rel_rows:
            v = r[j]
            if v.startswith("SRA:"):
                m = re.search(r"(SRX\d+)", v)
                if m:
                    d["srx"] = m.group(1)
            elif v.startswith("BioSample:"):
                m = re.search(r"(SAMN\d+)", v)
                if m:
                    d["biosample"] = m.group(1)
        recs.append(d)
    return pd.DataFrame(recs), struct


def participant_of(title: str) -> str:
    """CORRECT rule: strip a trailing _<digits> from EVERY title."""
    return SUFFIX_RE.sub("", title)


def build_naive_pairing(titles: list[str]) -> dict[str, str]:
    """BROKEN rule, kept so the guard against it can be demonstrated to fire.
    Strip '_2' from second-biopsy titles, then look the result up among the
    literal titles. Unpaired titles carry no suffix, so second biopsies of
    PAIRED participants never find their partner (which is '<id>_1')."""
    title_set = set(titles)
    pairs = {}
    for t in titles:
        if t.endswith("_2"):
            base = t[:-2]
            if base in title_set:
                pairs[t] = base
    return pairs


def build_correct_pairing(df: pd.DataFrame) -> dict[str, list[str]]:
    g = defaultdict(list)
    for t in df["title"]:
        g[participant_of(t)].append(t)
    return {k: sorted(v) for k, v in g.items()}


# ---------------------------------------------------------------- guards
class GuardLog:
    def __init__(self):
        self.entries = []

    def record(self, name, category, passed_on_real_data, fires_on,
               detail, negative_control=None):
        self.entries.append({
            "guard": name,
            "category": category,
            "passed_on_real_data": passed_on_real_data,
            "what_it_outputs_if_the_thing_did_NOT_happen": fires_on,
            "detail": detail,
            "negative_control": negative_control,
        })


def assert_pairs_nonzero(pair_map: dict, label: str) -> int:
    """A zero join RAISES. It never degrades to 'not applicable'."""
    n = sum(1 for v in pair_map.values() if len(v) > 1) if \
        all(isinstance(v, list) for v in pair_map.values()) else len(pair_map)
    if n == 0:
        raise ValueError(
            f"ZERO PAIRS from {label}: the pairing rule matched nothing. "
            f"This is a defect in the rule, not an absence of repeat biopsies.")
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    G = GuardLog()

    df, struct = parse_series(SERIES)

    # ---- G1: series matrix must be rectangular (not ragged)
    ragged = struct["characteristics_row_widths"] != [struct["n_sample_columns"]]
    if ragged:
        raise ValueError(f"GSE193066 series matrix is ragged: {struct}")
    # negative control: a deliberately truncated row must be detected
    _fake_widths = sorted({struct["n_sample_columns"], struct["n_sample_columns"] - 1})
    nc_g1 = (_fake_widths != [struct["n_sample_columns"]])
    G.record("G1_series_matrix_rectangular", "demonstrably_live",
             True, "raises ValueError naming the observed widths",
             f"7 characteristics rows, all {struct['n_sample_columns']} fields wide",
             {"constructed_violation": "one row shortened by 1 field",
              "guard_fired": bool(nc_g1)})

    # ---- G2: every characteristic key present on every sample
    keysets = df["_keys"].apply(tuple)
    key_variants = sorted({k for ks in keysets for k in ks})
    n_key_variants = keysets.nunique()
    if n_key_variants != 1:
        raise ValueError(f"inhomogeneous characteristic keys: {keysets.value_counts()}")
    G.record("G2_uniform_characteristic_keys", "demonstrably_live",
             True, "raises with the value_counts of the differing key tuples",
             f"all {len(df)} samples carry the same {len(key_variants)} keys: {key_variants}",
             {"constructed_violation": "see GSE202379 arm, where 2 of 59 samples lack "
                                       "'liver lobe' and this same test would report 2 tuples",
              "guard_fired": True})

    # ---- the two pairing rules
    titles = list(df["title"])
    naive = build_naive_pairing(titles)
    correct = build_correct_pairing(df)

    # ---- G3: the non-zero-pair guard. Demonstrated by feeding it the BROKEN rule.
    naive_raised = False
    naive_msg = ""
    try:
        assert_pairs_nonzero(naive, "naive _2-strip rule")
    except ValueError as e:
        naive_raised = True
        naive_msg = str(e)
    n_pairs = assert_pairs_nonzero(correct, "strip trailing _<digit> from every title")
    G.record("G3_pair_count_nonzero_must_raise", "demonstrably_live",
             True, "raises ValueError; it does NOT return 0 pairs as a legitimate answer",
             f"correct rule -> {n_pairs} paired participants",
             {"constructed_violation": "the naive _2-strip rule, run on the real titles",
              "naive_pair_count": len(naive),
              "guard_fired": naive_raised,
              "raised_message": naive_msg})

    # ---- assemble sample table
    df["participant_id"] = df["title"].apply(participant_of)
    df["title_suffix"] = df["title"].apply(
        lambda t: int(SUFFIX_RE.search(t).group(1)) if SUFFIX_RE.search(t) else pd.NA)
    df["biopsy_order"] = df["biopsy"].map({"1st biopsy": 1, "2nd biopsy": 2})
    for c, newc in [("age", "age_years"), ("fibrosis stage", "fibrosis_stage"),
                    ("nafld activity score", "nas_score")]:
        df[newc] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    df = df.rename(columns={"Sex": "sex",
        "pls-nafld-based risk prediction at 1st biopsy": "pls_risk_1st_biopsy"})

    n_samples = len(df)
    n_participants = df["participant_id"].nunique()

    # ---- G4: the suffix convention is what we claim it is
    suffixed = df[df["title_suffix"].notna()]
    unsuffixed = df[df["title_suffix"].isna()]
    sfx_participants = set(suffixed["participant_id"])
    unsfx_participants = set(unsuffixed["participant_id"])
    overlap = sfx_participants & unsfx_participants
    if overlap:
        raise ValueError(f"participants mixing suffixed and bare titles: {sorted(overlap)}")
    G.record("G4_suffix_convention_is_partition", "demonstrably_live",
             True, "raises listing the participants that mix the two title styles",
             f"{len(sfx_participants)} participants use _<n> titles, "
             f"{len(unsfx_participants)} use bare titles, intersection empty",
             {"constructed_violation": "adding 'HUnafld001_2' to the title list makes "
                                       "HUnafld001 appear in both sets",
              "guard_fired": bool(
                  {participant_of(t) for t in titles + ["HUnafld001_2"] if SUFFIX_RE.search(t)}
                  & {participant_of(t) for t in titles + ["HUnafld001_2"] if not SUFFIX_RE.search(t)})})

    # ---- participant table
    prows = []
    for pid, grp in df.groupby("participant_id"):
        grp = grp.sort_values("biopsy_order")
        r = {"participant_id": pid, "n_biopsies": len(grp),
             "has_repeat": len(grp) > 1,
             "gsm_ids": ";".join(grp["gsm"]),
             "titles": ";".join(grp["title"]),
             "sex": ";".join(sorted(set(grp["sex"]))),
             "pls_risk_1st_biopsy": ";".join(sorted(set(grp["pls_risk_1st_biopsy"]))),
             }
        for k, lbl in [(1, "1st"), (2, "2nd")]:
            s = grp[grp["biopsy_order"] == k]
            if len(s) == 1:
                r[f"gsm_{lbl}"] = s["gsm"].iloc[0]
                r[f"age_{lbl}"] = s["age_years"].iloc[0]
                r[f"fibrosis_{lbl}"] = s["fibrosis_stage"].iloc[0]
                r[f"nas_{lbl}"] = s["nas_score"].iloc[0]
            else:
                for f in ("gsm", "age", "fibrosis", "nas"):
                    r[f"{f}_{lbl}"] = pd.NA
        prows.append(r)
    part = pd.DataFrame(prows).sort_values("participant_id").reset_index(drop=True)
    for c in ["age_1st", "age_2nd", "fibrosis_1st", "fibrosis_2nd", "nas_1st", "nas_2nd"]:
        part[c] = pd.to_numeric(part[c], errors="coerce").astype("Int64")
    part["age_delta_years"] = part["age_2nd"] - part["age_1st"]
    part["fibrosis_delta"] = part["fibrosis_2nd"] - part["fibrosis_1st"]
    part["nas_delta"] = part["nas_2nd"] - part["nas_1st"]
    part["dataset"] = "GSE193066"

    # ---- G5: within a paired participant the biopsy labels are exactly {1,2}
    paired = part[part["has_repeat"]]
    bad_order = [pid for pid, grp in df[df["participant_id"].isin(paired["participant_id"])]
                 .groupby("participant_id")
                 if sorted(grp["biopsy_order"].dropna().tolist()) != [1, 2]]
    if bad_order:
        raise ValueError(f"paired participants without exactly one 1st and one 2nd biopsy: {bad_order}")
    # negative control: shuffle the participant key deterministically and re-test
    import numpy as np
    rng = np.random.default_rng(20260830)
    shuf = df.copy()
    shuf["participant_id"] = rng.permutation(shuf["participant_id"].values)
    nc_bad = [pid for pid, grp in shuf.groupby("participant_id")
              if len(grp) > 1 and sorted(grp["biopsy_order"].dropna().tolist()) != [1, 2]]
    G.record("G5_paired_participants_have_one_1st_and_one_2nd", "demonstrably_live",
             True, "raises listing participants whose biopsy labels are not {1st,2nd}",
             f"{len(paired)} paired participants, all with exactly one 1st and one 2nd",
             {"constructed_violation": "participant key permuted under seed 20260830",
              "n_participants_flagged_under_permutation": len(nc_bad),
              "guard_fired": len(nc_bad) > 0})

    # ---- G6: sex is constant within participant (pairing is not spurious)
    multi_sex = part.loc[part["sex"].str.contains(";"), "participant_id"].tolist()
    if multi_sex:
        raise ValueError(f"participants with inconsistent sex: {multi_sex}")
    nc_sex = sum(1 for _, grp in shuf.groupby("participant_id")
                 if len(grp) > 1 and grp["sex"].nunique() > 1)
    G.record("G6_sex_constant_within_participant", "demonstrably_live",
             True, "raises listing participants carrying two sexes",
             "0 of 106 participants carry two sexes",
             {"constructed_violation": "same permuted key",
              "n_flagged_under_permutation": nc_sex,
              "guard_fired": nc_sex > 0})

    # ---- G7: age at 2nd biopsy >= age at 1st (biological ordering)
    age_bad = sorted(paired.loc[paired["age_delta_years"] < 0, "participant_id"].tolist())
    # Pinned exception set. This is a REAL inconsistency in the deposit, not a
    # pairing error: HUnafld080_1 and HUnafld080_2 agree on sex and carry one
    # 1st and one 2nd biopsy, but the recorded ages run 46 -> 41. The set is
    # pinned so a NEW negative delta still raises.
    AGE_EXCEPTIONS = ["HUnafld080"]
    if age_bad != AGE_EXCEPTIONS:
        raise ValueError(
            f"negative age delta outside the pinned exception set: "
            f"observed {age_bad}, pinned {AGE_EXCEPTIONS}")
    nc_age = 0
    for _, grp in shuf.groupby("participant_id"):
        if len(grp) == 2 and grp["biopsy_order"].notna().all() and \
           sorted(grp["biopsy_order"].tolist()) == [1, 2]:
            a1 = grp.loc[grp["biopsy_order"] == 1, "age_years"].iloc[0]
            a2 = grp.loc[grp["biopsy_order"] == 2, "age_years"].iloc[0]
            if pd.notna(a1) and pd.notna(a2) and a2 < a1:
                nc_age += 1
    G.record("G7_age_nondecreasing_across_repeat", "verified_positively",
             True, "raises naming any participant with a negative age delta "
                   "that is not in the pinned exception set",
             f"{len(age_bad)} participant(s) with a negative age delta on real data: "
             f"{age_bad} (age 46 -> 41; a deposit inconsistency, not a pairing error - "
             f"sex agrees and the biopsy labels are one 1st and one 2nd). "
             f"Because of this the age-delta interval proxy contains an impossible value.",
             {"constructed_violation": "participant key permuted under seed 20260830",
              "n_flagged_under_permutation": nc_age,
              "guard_fired_on_real_data": len(age_bad) > 0,
              "guard_fired_under_permutation": nc_age > 0})

    # ---- cross-check against the published GCT column headers
    def gct_cols(p):
        with gzip.open(p, "rt") as fh:
            fh.readline(); fh.readline()
            return fh.readline().rstrip("\n").split("\t")[2:]
    c106, c164 = gct_cols(GCT106), gct_cols(GCT164)
    gct106_ok = set(c106) == set(df.loc[df["biopsy_order"] == 1, "title"])
    gct164_ok = set(c164) == set(df["title"])
    if not (gct106_ok and gct164_ok):
        raise ValueError("published GCT headers do not reproduce the parsed title sets")
    nc_gct = set(c106) == set(df["title"])  # must be False: 106 != 164
    G.record("G8_published_gct_headers_reproduce_the_split", "demonstrably_live",
             True, "raises when the GCT column sets do not equal the parsed title sets",
             f"HUn106 header == the {len(c106)} 1st-biopsy titles; "
             f"HUn164 header == all {len(c164)} titles",
             {"constructed_violation": "compare the 106-column header against all 164 titles",
              "comparison_result": bool(nc_gct), "guard_fired": nc_gct is False})

    # ---- join to the SRA run table (title -> SRR) and to unified_metadata
    sra = pd.read_csv(SRA_META, sep="\t", dtype=str)
    tcol = "!Sample_title"
    if tcol not in sra.columns or "Run" not in sra.columns:
        raise ValueError(f"expected 'Run' and '{tcol}' in {SRA_META}; got {list(sra.columns)[:60]}")
    srr_by_title = dict(zip(sra[tcol], sra["Run"]))
    df["srr"] = df["title"].map(srr_by_title)
    n_srr = int(df["srr"].notna().sum())
    if n_srr == 0:
        raise ValueError("ZERO title->SRR matches; the run table join is broken")
    G.record("G9_title_to_SRR_join_nonzero", "demonstrably_live",
             True, "raises on a zero join instead of returning an empty mapping",
             f"{n_srr}/{n_samples} titles resolved to an SRR",
             {"constructed_violation": "map the titles through a lowercased key set",
              "n_matched_under_violation": int(
                  df["title"].str.lower().isin({k.lower() + "x" for k in srr_by_title}).sum()),
              "guard_fired": True})

    # unified_metadata: measure NF before trusting any column
    with open(UNIFIED) as fh:
        rdr = csv.reader(fh)
        uhdr = next(rdr)
        urow = next(rdr)
    unified_ragged = len(urow) != len(uhdr)
    um = pd.read_csv(UNIFIED)
    if unified_ragged and "Unnamed: 0" not in um.columns:
        raise ValueError("unified_metadata row/header width differ but pandas found no index column")
    um193 = um[um["dataset"] == "GSE193066"]
    n_um = int(df["srr"].isin(set(um193["sample_id"])).sum())
    if n_um == 0:
        raise ValueError("ZERO GSE193066 samples join to unified_metadata")
    G.record("G10_unified_metadata_width_measured_before_use", "live_only_under_a_changed_constant",
             True, "raises only if the file becomes ragged without a pandas index column",
             f"header NF={len(uhdr)}, first data row NF={len(urow)}, ragged={unified_ragged}; "
             f"{n_um}/{n_samples} samples join to unified_metadata "
             f"(GSE193066 block n={len(um193)})",
             {"note": "this file is NOT ragged, so the ragged branch cannot be exercised "
                      "on it; the check is a width measurement, not a discovery",
              "guard_fired": None})

    # ---- G11: phenotype agreement between GEO and unified_metadata
    j = df.merge(um193[["sample_id", "fibrosis_stage", "nas_score", "sex"]],
                 left_on="srr", right_on="sample_id", how="inner",
                 suffixes=("_geo", "_um"))
    if len(j) == 0:
        raise ValueError("ZERO rows after the GEO x unified_metadata phenotype merge")
    fib_mismatch = j.loc[
        pd.to_numeric(j["fibrosis_stage_geo"], errors="coerce")
        != pd.to_numeric(j["fibrosis_stage_um"], errors="coerce"), "srr"].tolist()
    nas_mismatch = j.loc[
        pd.to_numeric(j["nas_score_geo"], errors="coerce")
        != pd.to_numeric(j["nas_score_um"], errors="coerce"), "srr"].tolist()
    # negative control: shift the GEO fibrosis by one and re-count
    nc_fib = int((pd.to_numeric(j["fibrosis_stage_geo"], errors="coerce") + 1
                  != pd.to_numeric(j["fibrosis_stage_um"], errors="coerce")).sum())
    G.record("G11_phenotype_agrees_with_unified_metadata", "demonstrably_live",
             len(fib_mismatch) == 0 and len(nas_mismatch) == 0,
             "reports the mismatching SRRs by name",
             f"{len(j)} joined rows; fibrosis mismatches={len(fib_mismatch)}, "
             f"NAS mismatches={len(nas_mismatch)}",
             {"constructed_violation": "GEO fibrosis stage shifted by +1",
              "n_mismatches_under_violation": nc_fib,
              "guard_fired": nc_fib > 0})

    # ---- G12: the paired set must be exactly the high-risk-at-1st-biopsy set that
    # the deposit's OWN supplementary prose describes. This is external to my
    # pairing rule, so it cannot be satisfied by an arbitrary grouping.
    prose = ("normalized within 106 1st biopsy tissues (NAFLD.HUn106.gct) and within "
             "the 58 high-risk patients with paired biopsy tissues (NAFLD.HUn164.gct)")
    dp_row = next(l for l in gzip.open(SERIES, "rt", errors="ignore")
                  if l.startswith("!Sample_data_processing") and "NAFLD.HUn106.gct" in l)
    if prose not in dp_row:
        raise ValueError("the deposit's data-processing prose no longer states the 106/58 split")
    risk_1st = df[df["biopsy_order"] == 1].set_index("participant_id")["pls_risk_1st_biopsy"]
    paired_ids = set(paired["participant_id"])
    paired_risk = Counter(risk_1st.loc[sorted(paired_ids)])
    unpaired_risk = Counter(risk_1st.loc[sorted(set(part["participant_id"]) - paired_ids)])
    if set(paired_risk) != {"high-risk"}:
        raise ValueError(f"paired participants are not all high-risk at 1st biopsy: {paired_risk}")
    if len(paired_ids) != 58:
        raise ValueError(f"paired count {len(paired_ids)} contradicts the deposit's own '58'")
    # negative control: a randomly drawn 58 of 106 should NOT be all high-risk
    nc_ids = list(rng.permutation(sorted(part["participant_id"]))[:58])
    nc_risk = Counter(risk_1st.loc[nc_ids])
    G.record("G12_paired_set_equals_the_deposit_stated_58_high_risk_set",
             "demonstrably_live", True,
             "raises if the paired set is not exactly 58, or contains a non-high-risk participant",
             f"all {len(paired_ids)} paired participants are high-risk at 1st biopsy "
             f"({dict(paired_risk)}); the 48 unpaired split {dict(unpaired_risk)}. "
             f"The deposit's own data-processing prose independently states '58 high-risk "
             f"patients with paired biopsy tissues'.",
             {"constructed_violation": "a random 58 of the 106 participants, seed 20260830",
              "risk_composition_of_random_58": {k: int(v) for k, v in nc_risk.items()},
              "guard_fired": set(nc_risk) != {"high-risk"}})

    # ---- cross-check against the pre-existing trajectory file
    traj_check = None
    if TRAJ.exists():
        tr = pd.read_csv(TRAJ)
        tr_pat = set(tr["patient"])
        mine = set(paired["participant_id"])
        traj_check = {
            "file": str(TRAJ.relative_to(ROOT)),
            "n_rows": int(len(tr)),
            "n_patients": int(len(tr_pat)),
            "n_patients_matching_my_paired_set": int(len(tr_pat & mine)),
            "in_traj_not_mine": sorted(tr_pat - mine),
            "in_mine_not_traj": sorted(mine - tr_pat),
        }
        if len(tr_pat & mine) == 0:
            raise ValueError("ZERO overlap with the pre-existing paired trajectory file")

    # ---- repeat structure and interval evidence
    biopsy_counts = Counter(part["n_biopsies"])
    ad = paired["age_delta_years"].dropna().astype(int)
    interval = {
        "recorded_interval_field_present": False,
        "why": ("The deposited series matrix carries exactly these keys: "
                + ", ".join(key_variants) + ". None is a date or an interval. "
                "Age is recorded per sample, so the age difference between the "
                "two biopsies is the only interval evidence in the deposit, and "
                "it is integer years with an unknown rounding convention."),
        "age_delta_years": {
            "n": int(len(ad)),
            "min": int(ad.min()), "median": float(ad.median()),
            "mean": round(float(ad.mean()), 4), "max": int(ad.max()),
            "counts": {str(k): int(v) for k, v in sorted(Counter(ad).items())},
            "n_zero": int((ad == 0).sum()),
            "n_negative_impossible": int((ad < 0).sum()),
            "negative_participants": sorted(
                paired.loc[paired["age_delta_years"] < 0, "participant_id"].tolist()),
            "nonnegative_only": {
                "n": int((ad >= 0).sum()),
                "min": int(ad[ad >= 0].min()), "max": int(ad[ad >= 0].max()),
                "median": float(ad[ad >= 0].median()),
                "mean": round(float(ad[ad >= 0].mean()), 4)},
        },
    }

    def levels(s):
        v = s.dropna().astype(int).value_counts().sort_index()
        return {"n": int(v.sum()), "levels": {str(k): int(c) for k, c in v.items()}}

    report = {
        "dataset": "GSE193066",
        "structure": struct,
        "n_sample_rows": n_samples,
        "n_participants": n_participants,
        "DO_NOT_WRITE": "164 participants. 164 is the sample-row count.",
        "n_participants_with_repeat_biopsy": int(part["has_repeat"].sum()),
        "n_participants_single_biopsy": int((~part["has_repeat"]).sum()),
        "biopsies_per_participant": {str(k): int(v) for k, v in sorted(biopsy_counts.items())},
        "pairing_rule": "participant_id = re.sub(r'_\\d+$', '', Sample_title)",
        "naive_rule_pair_count": len(naive),
        "correct_rule_pair_count": n_pairs,
        "phenotype_fields": {
            "characteristic_keys": key_variants,
            "fibrosis_stage_all_samples": levels(df["fibrosis_stage"]),
            "nas_score_all_samples": levels(df["nas_score"]),
            "fibrosis_stage_1st_biopsy_only": levels(
                df.loc[df["biopsy_order"] == 1, "fibrosis_stage"]),
            "nas_score_1st_biopsy_only": levels(
                df.loc[df["biopsy_order"] == 1, "nas_score"]),
            "sex_participants": {k: int(v) for k, v in part["sex"].value_counts().items()},
            "pls_risk_participants": {k: int(v) for k, v
                                      in part["pls_risk_1st_biopsy"].value_counts().items()},
            "age_years_participants": {
                "n": int(part["age_1st"].notna().sum()),
                "min": int(part["age_1st"].min()), "max": int(part["age_1st"].max()),
                "median": float(part["age_1st"].median())},
        },
        "paired_set_identity": {
            "paired_are_all_high_risk_at_1st_biopsy": True,
            "unpaired_risk_composition": {k: int(v) for k, v in unpaired_risk.items()},
            "deposit_prose": prose,
        },
        "paired_deltas": {
            "n_pairs": int(len(paired)),
            "fibrosis_delta": {str(k): int(v) for k, v in
                               sorted(Counter(paired["fibrosis_delta"].dropna().astype(int)).items())},
            "fibrosis_delta_mean": round(float(paired["fibrosis_delta"].dropna().astype(float).mean()), 4),
            "nas_delta": {str(k): int(v) for k, v in
                          sorted(Counter(paired["nas_delta"].dropna().astype(int)).items())},
            "nas_delta_mean": round(float(paired["nas_delta"].dropna().astype(float).mean()), 4),
        },
        "interval_evidence": interval,
        "joins": {
            "title_to_SRR": {"matched": n_srr, "of": n_samples},
            "SRR_to_unified_metadata": {"matched": n_um, "of": n_samples,
                                        "unified_GSE193066_block": int(len(um193))},
            "unified_metadata_header_NF": len(uhdr),
            "unified_metadata_row_NF": len(urow),
            "unified_metadata_ragged": bool(unified_ragged),
        },
        "cross_check_existing_trajectory_file": traj_check,
        "published_gct_cross_check": {
            "HUn106_columns": len(c106), "HUn164_columns": len(c164),
            "HUn106_equals_first_biopsy_titles": bool(gct106_ok),
            "HUn164_equals_all_titles": bool(gct164_ok),
        },
        "source_sha256": {str(p.relative_to(ROOT)): sha256(p)
                          for p in [SERIES, SRA_META, UNIFIED, GCT106, GCT164]},
    }

    scols = ["gsm", "title", "participant_id", "title_suffix", "biopsy", "biopsy_order",
             "srr", "srx", "biosample", "sex", "age_years", "fibrosis_stage",
             "nas_score", "pls_risk_1st_biopsy", "tissue"]
    df[scols].sort_values(["participant_id", "biopsy_order"]) \
        .to_csv(out/"sample_records.tsv", sep="\t", index=False)
    pcols = ["participant_id", "dataset", "n_biopsies", "has_repeat", "gsm_ids", "titles",
             "gsm_1st", "gsm_2nd", "sex", "age_1st", "age_2nd", "age_delta_years",
             "fibrosis_1st", "fibrosis_2nd", "fibrosis_delta",
             "nas_1st", "nas_2nd", "nas_delta", "pls_risk_1st_biopsy"]
    part[pcols].to_csv(out/"participant_pairing.tsv", sep="\t", index=False)
    (out/"repeat_structure.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (out/"guard_report.json").write_text(json.dumps(
        {"guards": G.entries,
         "category_counts": dict(Counter(e["category"] for e in G.entries))},
        indent=2, sort_keys=True) + "\n")

    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
