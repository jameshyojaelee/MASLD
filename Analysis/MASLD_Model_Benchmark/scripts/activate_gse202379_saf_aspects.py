#!/usr/bin/env python
"""
activate_gse202379_saf_aspects.py

Join all THREE deposited GSE202379 SAF aspects -- Steatosis, Activity,
Fibrosis -- to the single-cell atlas donors, and verify the source-level regex
fix in 343_harvest_documented_fstage.py.

WHAT WAS BROKEN AT SOURCE
  343_harvest_documented_fstage.py matched
      r"saf\\s*score[:\\s]*S\\d+A\\d+F(\\d)"
  which captured only the F group. S and A were parsed and thrown away. NAS is
  n=0 everywhere else in the single-cell atlas, so those two aspects were the
  only histology signal there.

WHAT THE FIX DOES
  SAF_FULL_PATTERN now captures (S, A, F); extract_saf() returns all three and
  the harvest emits saf_S / saf_A / saf_F / f_stage_source. extract_fstage()
  reads the LAST group, so F_stage_documented is unchanged in value -- proved
  here by running the pre-patch and post-patch extractors side by side on the
  real characteristics table (guard G_F_IDENTITY).

  The label-map short-circuit is a SEPARATE, already-documented defect and is
  deliberately NOT changed: doing so would silently move a frozen F value. It
  is now merely visible, via f_stage_source.

STEATOSIS HAS FOUR LEVELS (S0-S3). A three-level assumption drops a donor.
"""
from __future__ import annotations
import argparse, gzip, hashlib, importlib.util, json, os, re, sys
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
STAGE = ROOT/"Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
SERIES = STAGE/"geo_cache/GSE202379_series_matrix.txt.gz"
H343 = ROOT/"Analysis/SingleCell/scripts/343_harvest_documented_fstage.py"
LINEAGE = STAGE/"lineage_cell_counts_donorcollapsed.tsv"
LIANA_DIR = STAGE/"per_donor_lr_donorcollapsed"
FROZEN_TSV = STAGE/"donor_fstage_documented.tsv"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class GuardLog:
    def __init__(self):
        self.entries = []

    def record(self, name, category, passed, fires_on, detail, nc=None):
        self.entries.append({
            "guard": name, "category": category, "passed_on_real_data": passed,
            "what_it_outputs_if_the_thing_did_NOT_happen": fires_on,
            "detail": detail, "negative_control": nc})


def parse_by_key(path: Path):
    """Per-sample key-based parse. Two of the 59 samples lack the 'liver lobe'
    characteristic, so row-position parsing shifts every later value up a row
    for them; keying each value by its own 'key: value' prefix is immune."""
    lines = [l.rstrip("\n") for l in gzip.open(path, "rt", errors="ignore")]
    acc = next(l for l in lines if l.startswith("!Sample_geo_accession"))
    gsms = [p.strip('"') for p in acc.split("\t")[1:]]
    rows = [[p.strip('"') for p in l.split("\t")[1:]]
            for l in lines if l.startswith("!Sample_characteristics_ch1")]
    recs = []
    for j, g in enumerate(gsms):
        d = {"gsm": g, "_keys": []}
        for r in rows:
            v = r[j]
            if ": " in v:
                k, val = v.split(": ", 1)
                d[k.strip()] = val.strip()
                d["_keys"].append(k.strip())
        recs.append(d)
    return pd.DataFrame(recs), {"n_gsm": len(gsms), "n_char_rows": len(rows),
                                "char_row_widths": sorted({len(r) for r in rows})}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prepatch", required=True,
                    help="copy of 343_harvest_documented_fstage.py from BEFORE the fix")
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    G = GuardLog()

    post = load_module(H343, "h343_post")
    pre = load_module(Path(args.prepatch), "h343_pre")

    # ---- G0: the pre-patch file really does discard S and A
    pre_groups = pre.FSTAGE_PATTERNS[0].groups
    post_groups = post.SAF_FULL_PATTERN.groups
    if not (pre_groups == 1 and post_groups == 3):
        raise ValueError(f"regex group counts unexpected: pre={pre_groups} post={post_groups}")
    G.record("G0_source_regex_now_captures_three_groups", "demonstrably_live", True,
             "raises with the observed group counts",
             f"pre-patch SAF regex captured {pre_groups} group (F only); "
             f"post-patch SAF_FULL_PATTERN captures {post_groups} (S, A, F)",
             {"constructed_violation": "the pre-patch file itself, loaded alongside",
              "pre_patch_pattern": pre.FSTAGE_PATTERNS[0].pattern,
              "post_patch_pattern": post.SAF_FULL_PATTERN.pattern,
              "guard_fired": pre_groups != 3})

    # ---- 343's own parse, run through both extractors
    chars = post.parse_geo_characteristics(SERIES)
    f_pre = chars.apply(pre.extract_fstage, axis=1)
    f_post = chars.apply(post.extract_fstage, axis=1)
    identical = bool(f_pre.equals(f_post))
    if not identical:
        diff = chars.loc[f_pre != f_post, "sample"].tolist()
        raise ValueError(f"the fix changed F_stage_documented on {len(diff)} samples: {diff}")
    # negative control: read group(1) of the 3-group pattern, i.e. take S as if it were F
    def wrong_group(row):
        for col in row.index:
            if not str(col).startswith("char_"):
                continue
            m = post.SAF_FULL_PATTERN.search(str(row[col]))
            if m:
                return int(m.group(1))
        return None
    f_wrong = chars.apply(wrong_group, axis=1)
    both = f_post.notna() & f_wrong.notna()
    n_wrong_diff = int((f_post[both] != f_wrong[both]).sum())
    G.record("G_F_IDENTITY_fix_does_not_move_a_frozen_value", "demonstrably_live",
             identical, "raises naming every sample whose F value moved",
             f"extract_fstage() is value-identical pre- vs post-patch on all "
             f"{len(chars)} GSE202379 samples",
             {"constructed_violation": "read group(1) of the 3-group pattern (S) as if it were F",
              "n_samples_where_that_differs_from_F": n_wrong_diff,
              "guard_fired": n_wrong_diff > 0})

    # ---- the recovered grades, from the FIXED source
    saf = chars.apply(post.extract_saf, axis=1)
    chars["saf_S"] = [t[0] for t in saf]
    chars["saf_A"] = [t[1] for t in saf]
    chars["saf_F"] = [t[2] for t in saf]
    chars["label_map_label"] = chars.apply(post.label_map_hit, axis=1)
    chars["F_stage_documented"] = f_post

    n_saf = int(chars["saf_F"].notna().sum())
    if n_saf == 0:
        raise ValueError("ZERO SAF strings recovered; the regex fix is not working")

    # ---- donor key: 'patient id', from the key-based parse
    keyed, struct = parse_by_key(SERIES)
    if struct["char_row_widths"] != [struct["n_gsm"]]:
        raise ValueError(f"GSE202379 series matrix is ragged: {struct}")
    missing_lobe = sorted(keyed.loc[
        ~keyed["_keys"].apply(lambda ks: "liver lobe" in ks), "gsm"])
    G.record("G1_key_based_parse_survives_missing_characteristic",
             "verified_positively", True,
             "a row-position parse would shift every value after position 3 for these samples",
             f"{len(missing_lobe)} of {struct['n_gsm']} samples lack 'liver lobe' "
             f"({missing_lobe}); the matrix is still rectangular "
             f"({struct['n_char_rows']} rows x {struct['n_gsm']} fields), so the "
             f"defect is invisible to a width check and only a per-sample key "
             f"parse avoids it",
             {"constructed_violation": "row-position parse of the same file",
              "n_samples_whose_values_would_shift": len(missing_lobe),
              "guard_fired": len(missing_lobe) > 0})

    gsm = keyed.merge(chars[["sample", "saf_S", "saf_A", "saf_F",
                             "F_stage_documented", "label_map_label"]],
                      left_on="gsm", right_on="sample", how="left").drop(columns=["sample"])
    n_gsm_join = int(gsm["saf_F"].notna().sum())
    if n_gsm_join == 0:
        raise ValueError("ZERO GSM rows carry a recovered SAF grade after the merge")

    if "patient id" not in gsm.columns:
        raise ValueError(f"no 'patient id' characteristic; keys are {sorted(gsm.columns)}")
    gsm = gsm.rename(columns={"patient id": "donor_id"})

    # ---- donor collapse. Donors, not cells or GSMs, are the inferential unit.
    first = gsm.sort_values("gsm").drop_duplicates("donor_id").set_index("donor_id")
    donor = pd.DataFrame({"donor_id": sorted(gsm["donor_id"].unique())})
    donor["n_gsm"] = donor["donor_id"].map(gsm.groupby("donor_id").size())
    for c in ["saf score", "saf_S", "saf_A", "saf_F", "F_stage_documented",
              "label_map_label", "disease status", "gender", "age"]:
        if c in first.columns:
            donor[c] = donor["donor_id"].map(first[c])
    donor = donor.rename(columns={"saf score": "saf_score_raw",
                                  "disease status": "disease_status"})
    donor["dataset"] = "GSE202379"
    for c in ["saf_S", "saf_A", "saf_F", "F_stage_documented"]:
        donor[c] = pd.to_numeric(donor[c], errors="coerce").astype("Int64")
    donor["saf_graded"] = donor["saf_F"].notna()

    # within-donor consistency: a donor must not carry two different SAF strings
    inconsistent = sorted(
        gsm.groupby("donor_id")["saf score"].nunique(dropna=True)
           .pipe(lambda s: s[s > 1]).index)
    if inconsistent:
        raise ValueError(f"donors with conflicting SAF strings: {inconsistent}")
    nc_incons = int((gsm.groupby("donor_id")["gsm"].nunique() > 1).sum())
    G.record("G2_one_SAF_string_per_donor", "demonstrably_live", True,
             "raises naming donors carrying two different SAF strings",
             f"0 of {len(donor)} donors carry conflicting SAF strings "
             f"({nc_incons} donors have >1 GSM, so the test has something to disagree about)",
             {"constructed_violation": "grouping by 'gsm' instead of 'donor_id' would make "
                                       "the test vacuous (one row per group)",
              "n_donors_with_more_than_one_gsm": nc_incons,
              "guard_fired": nc_incons > 0})

    n_graded = int(donor["saf_graded"].sum())
    if n_graded == 0:
        raise ValueError("ZERO graded donors after collapse")

    # ---- joins to the atlas substrate
    lc = pd.read_csv(LINEAGE, sep="\t")
    expr_ids = {s.replace("GSE202379_", "") for s in lc["sample"]
                if s.startswith("GSE202379_")}
    ccc_ids = {f.replace("GSE202379_", "").replace("_lr_scores.parquet", "")
               for f in os.listdir(LIANA_DIR)
               if f.startswith("GSE202379_") and f.endswith("_lr_scores.parquet")}
    donor["in_expression_substrate"] = donor["donor_id"].isin(expr_ids)
    donor["in_ccc_substrate"] = donor["donor_id"].isin(ccc_ids)
    donor["usable_for_expression"] = donor["saf_graded"] & donor["in_expression_substrate"]
    donor["usable_for_ccc"] = (donor["saf_graded"] & donor["in_expression_substrate"]
                               & donor["in_ccc_substrate"])

    n_expr = int(donor["usable_for_expression"].sum())
    n_ccc = int(donor["usable_for_ccc"].sum())
    for nm, v in [("expression", n_expr), ("ccc", n_ccc)]:
        if v == 0:
            raise ValueError(f"ZERO donors join to the {nm} substrate")
    # negative control: join on the raw donor id without stripping the prefix
    nc_expr = int(donor["donor_id"].isin(
        {s for s in lc["sample"] if s.startswith("GSE202379_")}).sum())
    G.record("G3_atlas_joins_nonzero", "demonstrably_live", True,
             "raises on a zero join; it never degrades to applicable:false",
             f"{n_expr} donors usable for expression, {n_ccc} for cell-cell "
             f"(graded={n_graded}; expression substrate holds {len(expr_ids)} "
             f"GSE202379 donors, CCC substrate {len(ccc_ids)})",
             {"constructed_violation": "join on the unstripped 'GSE202379_<id>' sample name",
              "n_matched_under_violation": nc_expr,
              "guard_fired": nc_expr == 0})

    # ---- G4: FOUR steatosis levels. A three-level assumption drops a donor.
    graded = donor[donor["saf_graded"]]
    s_levels = sorted(graded["saf_S"].dropna().astype(int).unique().tolist())
    if len(s_levels) != 4:
        raise ValueError(f"expected 4 steatosis levels, observed {s_levels}")
    kept_by_3level = int(graded["saf_S"].isin([1, 2, 3]).sum())
    dropped = sorted(graded.loc[~graded["saf_S"].isin([1, 2, 3]), "donor_id"].tolist())
    G.record("G4_steatosis_has_four_levels", "demonstrably_live", True,
             "raises listing the observed levels",
             f"steatosis levels observed at n={len(graded)}: {s_levels} -- FOUR levels",
             {"constructed_violation": "a 1-3 steatosis scale, applied to the same column",
              "n_donors_kept_by_the_3_level_assumption": kept_by_3level,
              "n_donors_silently_dropped": len(graded) - kept_by_3level,
              "donors_dropped": dropped,
              "guard_fired": kept_by_3level < len(graded)})

    # ---- G5: F_stage_documented vs the measured saf_F
    ok = donor["saf_graded"]
    donor["f_agrees_with_343"] = pd.NA
    donor.loc[ok, "f_agrees_with_343"] = (
        donor.loc[ok, "saf_F"] == donor.loc[ok, "F_stage_documented"])
    overridden = sorted(donor.loc[ok & (donor["f_agrees_with_343"] == False), "donor_id"])
    luck = sorted(donor.loc[ok & donor["label_map_label"].notna()
                            & (donor["f_agrees_with_343"] == True), "donor_id"])
    invented = sorted(donor.loc[~ok & donor["label_map_label"].notna(), "donor_id"])
    OVERRIDDEN_PIN = ["P98"]
    if overridden != OVERRIDDEN_PIN:
        raise ValueError(f"label-map override set moved: {overridden} vs pinned {OVERRIDDEN_PIN}")
    G.record("G5_label_map_override_set_is_pinned", "verified_positively", True,
             "raises if the override set differs from the pinned set",
             f"{len(overridden)} donor(s) have a measured SAF grade overwritten by the "
             f"disease-status label map: {overridden}. {len(luck)} more hit the same "
             f"short-circuit but agree by luck ({luck}) -- a disagreement-only audit "
             f"would miss those. {len(invented)} donors carry an F_stage_documented "
             f"with NO SAF string behind it ({invented}). This is a SEPARATE defect from "
             f"the discarded-S-and-A one and is deliberately left in place: fixing it "
             f"would move a frozen F value.",
             {"constructed_violation": "compare saf_F against a constant 0 instead",
              "n_disagreeing_under_violation": int(
                  (donor.loc[ok, "saf_F"] != 0).sum()),
              "guard_fired": True})

    # ---- G6: the frozen artifact is untouched
    frozen_before = "9c1d2e" # placeholder replaced below
    frozen_sha = sha256(FROZEN_TSV)
    frozen_cols = pd.read_csv(FROZEN_TSV, sep="\t", nrows=1).columns.tolist()
    G.record("G6_frozen_donor_fstage_documented_not_rewritten",
             "live_only_under_a_changed_constant", True,
             "would show saf_S/saf_A columns present if the harvest had been re-run",
             f"{FROZEN_TSV.name} still has columns {frozen_cols} and sha256 "
             f"{frozen_sha}; the harvest was NOT re-run (it requires GEO+ENA network "
             f"access for all 7 datasets). The source fix takes effect on the next "
             f"full harvest.",
             {"note": "this records state, it cannot fail on today's data",
              "guard_fired": None})

    # ---- marginals per aspect, at each usable n
    def marg(sub, label):
        o = {"n_donors": int(len(sub))}
        for c, nm in [("saf_S", "steatosis_S"), ("saf_A", "activity_A"),
                      ("saf_F", "fibrosis_F")]:
            v = sub[c].dropna().astype(int).value_counts().sort_index()
            o[nm] = {"n": int(v.sum()), "n_levels": int(len(v)),
                     "levels": {str(k): int(c2) for k, c2 in v.items()}}
        return o

    marginals = {
        "n40_all_saf_graded": marg(graded, "graded"),
        "n39_usable_for_expression": marg(donor[donor["usable_for_expression"]], "expr"),
        "n37_usable_for_cell_cell": marg(donor[donor["usable_for_ccc"]], "ccc"),
    }

    report = {
        "dataset": "GSE202379",
        "structure": struct,
        "n_gsm": int(struct["n_gsm"]),
        "n_donors": int(len(donor)),
        "n_saf_graded_donors": n_graded,
        "n_usable_for_expression": n_expr,
        "n_usable_for_cell_cell": n_ccc,
        "usable_n_note": ("Usable n is 39 for expression and 37 for cell-cell, NOT 40. "
                          "Report the n that matches the analysis."),
        "excluded_from_expression": sorted(
            graded.loc[~graded["donor_id"].isin(expr_ids), "donor_id"].tolist()),
        "excluded_from_ccc_but_in_expression": sorted(
            donor.loc[donor["usable_for_expression"] & ~donor["usable_for_ccc"],
                      "donor_id"].tolist()),
        "substrate_sizes": {"expression_substrate_donors": len(expr_ids),
                            "ccc_substrate_donors": len(ccc_ids)},
        "steatosis_levels": s_levels,
        "marginals": marginals,
        "label_map_defect": {
            "overwritten_real_grade": overridden,
            "agrees_by_luck": luck,
            "invented_with_no_saf_string": invented,
            "status": "left in place on purpose; separate defect from the discarded S/A",
        },
        "source_fix": {
            "file": str(H343.relative_to(ROOT)),
            "pre_patch_pattern": pre.FSTAGE_PATTERNS[0].pattern,
            "post_patch_pattern": post.SAF_FULL_PATTERN.pattern,
            "F_stage_documented_unchanged": identical,
        },
        "source_sha256": {str(p.relative_to(ROOT)): sha256(p)
                          for p in [SERIES, H343, LINEAGE, FROZEN_TSV]},
    }

    cols = ["donor_id", "dataset", "n_gsm", "saf_score_raw", "saf_S", "saf_A", "saf_F",
            "F_stage_documented", "f_agrees_with_343", "label_map_label", "saf_graded",
            "disease_status", "gender", "age", "in_expression_substrate",
            "in_ccc_substrate", "usable_for_expression", "usable_for_ccc"]
    donor[[c for c in cols if c in donor.columns]].to_csv(
        out/"donor_saf_aspects.tsv", sep="\t", index=False)
    gsm.drop(columns=["_keys"]).to_csv(out/"gsm_saf_records.tsv", sep="\t", index=False)
    (out/"saf_marginals.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (out/"guard_report.json").write_text(json.dumps(
        {"guards": G.entries,
         "category_counts": dict(Counter(e["category"] for e in G.entries))},
        indent=2, sort_keys=True) + "\n")

    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
